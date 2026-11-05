"""Generate `lbs_reference.h5`: pointings, HWP angles, hit maps and h-maps
computed with litebird_sim, against which Falcons is tested.

Run from the repository root with a Python environment that has
litebird_sim installed (numpy >= 2):

    python test/reference/generate_lbs_reference.py

All inputs are synthetic. The litebird_sim commit and the versions of the
numerical libraries are stored as attributes of the output file.

Layout (h5py, row-major; Julia reads every array with its axes reversed):

    config/                       scan, instrument and detector parameters
    <case>/pointings_ecl          (ndet, nsamp_stored, 3)   θ, φ, ψ in rad
    <case>/pointings_gal          (ndet, nsamp_stored, 3)
    <case>/hwp_angle              (nsamp_stored,)
    <case>/sample_index           (nsamp_stored,)           0-based sample index
    E_gal/, E_ecl/                h-maps: <det>/<n>,<m>/{Re,Im}
    F/hitmap                      (npix,)
    G/, H/, I/                    time-dependent inputs stored next to the outputs
    imo/                          values litebird_sim reads from mock_imo/schema.json
    J/                            pointings of an observation built from the mock IMO

The script also (re)writes the synthetic instrument database
`mock_imo/schema.json` used by the IMO tests.
"""

import json
import os
import subprocess
import tempfile
import uuid

import astropy
import ducc0
import h5py
import numpy as np

import litebird_sim as lbs
from litebird_sim.coordinates import ECL_TO_GAL_QUAT, ECL_TO_GAL_ROT_MATRIX

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "lbs_reference.h5")
MOCK_IMO_DIR = os.path.join(HERE, "mock_imo")

SPIN_SUN_ANGLE_RAD = np.deg2rad(45.0)
PRECESSION_RATE_HZ = 1.0 / (192.348 * 60.0)
SPIN_RATE_HZ = 0.05 / 60.0
SPIN_BORESIGHT_ANGLE_RAD = np.deg2rad(50.0)
SAMPLING_RATE_HZ = 19.0
DURATION_S = 600.0
START_TIME_B = 1.234e6
HWP_ANG_SPEED_RADPSEC = 46.0 * 2.0 * np.pi / 60.0
HWP_START_ANGLE_RAD = 0.2
STRIDE = 10  # store every STRIDE-th sample of the pointings
NSIDE = 32
N_M_COUPLES = np.array([[0, 0], [1, 0], [2, 0], [4, 0], [0, 4], [2, 4], [2, -4]])

# (theta, phi, psi) of the detectors of case B, in degrees
DETECTOR_ANGLES_DEG = [(0.0, 0.0, 0.0), (5.0, 30.0, 45.0), (2.0, 200.0, -20.0)]
POL_ANGLES_RAD = [0.0, 0.5, 1.1]


def det_quat(theta, phi, psi=0.0):
    """Rz(phi) · Ry(theta) · Rz(psi), (x, y, z, w)."""
    q = np.array(lbs.quat_rotation_z(psi))
    lbs.quat_left_multiply(q, *lbs.quat_rotation_y(theta))
    lbs.quat_left_multiply(q, *lbs.quat_rotation_z(phi))
    return q


def scanning_strategy():
    return lbs.SpinningScanningStrategy(
        spin_sun_angle_rad=SPIN_SUN_ANGLE_RAD,
        precession_rate_hz=PRECESSION_RATE_HZ,
        spin_rate_hz=SPIN_RATE_HZ,
    )


def instrument_a():
    return lbs.InstrumentInfo(name="instr_a", spin_boresight_angle_rad=SPIN_BORESIGHT_ANGLE_RAD)


def instrument_b():
    return lbs.InstrumentInfo(
        name="instr_b",
        spin_boresight_angle_rad=SPIN_BORESIGHT_ANGLE_RAD,
        spin_rotangle_rad=np.pi,
        boresight_rotangle_rad=0.3,
    )


def detectors_b():
    return [
        lbs.DetectorInfo(
            name=f"det{i}",
            sampling_rate_hz=SAMPLING_RATE_HZ,
            pol_angle_rad=POL_ANGLES_RAD[i],
            quat=det_quat(*np.deg2rad(angles)),
        )
        for i, angles in enumerate(DETECTOR_ANGLES_DEG)
    ]


def hwp_b():
    return lbs.IdealHWP(ang_speed_radpsec=HWP_ANG_SPEED_RADPSEC, start_angle_rad=HWP_START_ANGLE_RAD)


def make_sim(workdir, name, start_time, instrument, detectors, delta_time_s, hwp=None,
             num_of_obs=1, setup_attitude=None):
    sim = lbs.Simulation(
        base_path=os.path.join(workdir, name),
        start_time=start_time,
        duration_s=DURATION_S,
        random_seed=12345,
        imo=lbs.Imo(flatfile_location=str(lbs.PTEP_IMO_LOCATION)),
    )
    sim.set_instrument(instrument)
    sim.create_observations(detectors=detectors, num_of_obs_per_detector=num_of_obs, allocate_tod=False)
    sim.set_scanning_strategy(scanning_strategy=scanning_strategy(), delta_time_s=delta_time_s)
    if setup_attitude is not None:
        setup_attitude(sim)
    if hwp is not None:
        sim.set_hwp(hwp)
    sim.prepare_pointings()
    return sim


def write(group, name, data):
    data = np.asarray(data)
    if data.ndim == 0:
        group.create_dataset(name, data=data)
    else:
        group.create_dataset(name, data=data, compression="gzip", compression_opts=9, shuffle=True)


def store_pointings(f, case, obs, with_hwp):
    g = f.create_group(case)
    idx = np.arange(0, obs.n_samples_process, STRIDE)
    ptg_ecl, hwp_angle = obs.get_pointings("all", compute_hwp_angle=with_hwp)
    ptg_gal, _ = obs.get_pointings(
        "all", compute_hwp_angle=False, output_coordinate_system=lbs.CoordinateSystem.Galactic
    )
    write(g, "pointings_ecl", ptg_ecl[:, idx, :])
    write(g, "pointings_gal", ptg_gal[:, idx, :])
    write(g, "sample_index", idx)
    if with_hwp:
        write(g, "hwp_angle", hwp_angle[idx])
    g.attrs["start_time_s"] = float(obs.start_time_process)
    g.attrs["n_samples"] = int(obs.n_samples_process)
    return g


def store_h_maps(f, case, sim, coordinate_system):
    res = lbs.make_h_maps(
        observations=sim.observations,
        nside=NSIDE,
        n_m_couples=N_M_COUPLES,
        output_coordinate_system=coordinate_system,
        save_to_file=False,
    )
    g = f.create_group(case)
    g.attrs["coordinate_system"] = coordinate_system.name
    g.attrs["nside"] = NSIDE
    write(g, "n_m_couples", N_M_COUPLES)
    for det, maps in res.h_maps.items():
        for (n, m), hmap in maps.items():
            write(g, f"{det}/{n},{m}/Re", hmap.real)
            write(g, f"{det}/{n},{m}/Im", hmap.imag)


def wobble_series(start_time, nquats, rate_hz, amplitude_rad, period_s, nominal):
    """Rx(a sin(2π t / P)) · nominal sampled at start_time + k / rate_hz."""
    t = start_time + np.arange(nquats) / rate_hz
    quats = np.empty((nquats, 4))
    for k in range(nquats):
        q = np.array(nominal, dtype=float)
        lbs.quat_left_multiply(q, *lbs.quat_rotation_x(amplitude_rad * np.sin(2 * np.pi * (t[k] - start_time) / period_s)))
        quats[k] = q
    return quats


# ----------------------------------------------------------------------------
# Synthetic instrument database (litebird_sim / libinsdb flat-file layout)


def _uuid(path):
    return str(uuid.uuid5(uuid.NAMESPACE_URL, "falcons-mock-imo" + path))


def _unit(q):
    q = np.asarray(q, dtype=float)
    return q / np.linalg.norm(q)


MOCK_RELEASES = ("vMOCK1", "vMOCK2")
MOCK_TIME_DEP_START = 99.0
MOCK_CHANNELS = {
    "M-100": ["M100_00_A", "M100_00_B", "M100_01_A"],
    "M-200": ["M200_00_A"],
}


def _detector_metadata(name, channel, k):
    angles = np.deg2rad([1.0 + k, 30.0 * k, 10.0 * k])
    meta = {
        "name": name,
        "wafer": f"W{k % 2}",
        "pixel": str(k + 7),
        "pixtype": "MP1",
        "channel": channel,
        "squid": k + 1,
        "sampling_rate_hz": 19,
        "fwhm_arcmin": 30.0 + k,
        "ellipticity": 1.05,
        "psi_rad": 0.1 * k,
        "bandcenter_ghz": float(channel[2:]),
        "bandwidth_ghz": 0.23 * float(channel[2:]),
        "net_ukrts": 50.0 + k,
        "pol_sensitivity_ukarcmin": 70.0 + k,
        "fknee_mhz": 20,
        "fmin_hz": 1e-5,
        "alpha": 1.0,
        "pol": "T" if k % 2 == 0 else "B",
        "orient": "Q",
        "quat": list(2.0 * det_quat(*angles)),  # not normalised on purpose
        "pol_angle_rad": 0.25 * k,
        "pointing_theta_phi_psi_deg": list(np.rad2deg(angles)),
        "pointing_u_v": [0.1 * k, -0.2 * k],
        "unrecognised_field": 1.0,  # not a DetectorInfo field: ignored
    }
    if name == "M100_01_A":
        t = MOCK_TIME_DEP_START + np.arange(70.0)
        quats = []
        for tk in t:
            q = det_quat(*angles)
            lbs.quat_left_multiply(q, *lbs.quat_rotation_x(np.deg2rad(2.0 / 60.0) * np.sin(2 * np.pi * (tk - t[0]) / 31.0)))
            quats.append(list(q))
        meta["quat"] = {"start_time": MOCK_TIME_DEP_START, "sampling_rate_hz": 1.0, "quats": quats}
    return meta


def build_mock_imo():
    entities, quantities, data_files = [], [], []
    release_files = {tag: [] for tag in MOCK_RELEASES}

    def add_quantity(entity_path, name, metadata, releases=MOCK_RELEASES, date="2026-01-01T00:00:00+00:00", key=""):
        q_uuid = _uuid(entity_path + "/" + name)
        if not any(q["uuid"] == q_uuid for q in quantities):
            quantities.append({"uuid": q_uuid, "name": name, "format_spec": _uuid("/format"), "entity": _uuid(entity_path)})
        f_uuid = _uuid(entity_path + "/" + name + "/" + key)
        data_files.append({"uuid": f_uuid, "name": name, "upload_date": date, "metadata": metadata,
                           "quantity": q_uuid, "spec_version": "1.0", "dependencies": []})
        for tag in releases:
            release_files[tag].append(f_uuid)
        return f_uuid

    channel_entities = []
    channel_uuids = []
    for channel, dets in MOCK_CHANNELS.items():
        det_entities = [{"uuid": _uuid(f"/satellite/MOCKFT/{channel}/{d}"), "name": d} for d in dets]
        channel_entities.append({"uuid": _uuid(f"/satellite/MOCKFT/{channel}"), "name": channel, "children": det_entities})
        det_uuids = []
        for k, d in enumerate(dets):
            path = f"/satellite/MOCKFT/{channel}/{d}"
            if d == "M100_00_A":
                old = _detector_metadata(d, channel, k)
                old["quat"] = [0.0, 0.0, 0.0, 1.0]
                add_quantity(path, "detector_info", old, releases=("vMOCK1",), date="2025-01-01T00:00:00+00:00", key="old")
                det_uuids.append(add_quantity(path, "detector_info", _detector_metadata(d, channel, k), releases=("vMOCK2",)))
            else:
                det_uuids.append(add_quantity(path, "detector_info", _detector_metadata(d, channel, k)))
        center = float(channel[2:])
        ch_meta = {
            "channel": channel, "bandcenter_ghz": center, "bandwidth_ghz": 0.23 * center,
            "sampling_rate_hz": 19, "fwhm_arcmin": 30.0, "fknee_mhz": 20, "fmin_hz": 1e-5, "alpha": 1.0,
            "detector_names": dets,
            # one entry given as a URL to exercise the UUID extraction
            "detector_objs": [f"/data_files/{det_uuids[0]}"] + det_uuids[1:],
            "pol_sensitivity_channel_ukarcmin": 20.0,
        }
        if channel == "M-100":
            ch_meta.update(number_of_detectors=len(dets), net_detector_ukrts=50.0, net_channel_ukrts=0.0)
        else:
            ch_meta.update(net_detector_ukrts=0.0, net_channel_ukrts=10.0)
        channel_uuids.append(add_quantity(f"/satellite/MOCKFT/{channel}", "channel_info", ch_meta))

    add_quantity("/satellite/MOCKFT", "instrument_info", {
        "name": "MOCKFT", "boresight_rotangle_deg": 12.0, "spin_boresight_angle_deg": 50.0,
        "spin_rotangle_deg": 180.0, "hwp_rpm": 46.0, "number_of_channels": len(MOCK_CHANNELS),
        "channel_names": list(MOCK_CHANNELS), "channel_objs": channel_uuids, "sampling_rate_hz": 19,
        "wafers": ["W0", "W1"], "waferspace": 8.0,
    })
    add_quantity("/satellite", "scanning_parameters", {
        "spin_sun_angle_deg": 45.0, "precession_period_min": 192.348, "spin_rate_rpm": 0.05,
        "mission_duration_year": 3.0,
    })
    entities = [{"uuid": _uuid("/satellite"), "name": "satellite", "children": [
        {"uuid": _uuid("/satellite/MOCKFT"), "name": "MOCKFT", "children": channel_entities}]}]
    schema = {
        "instrumentdb": {"version": "mock"},
        "entities": entities,
        "format_specifications": [{"uuid": _uuid("/format"), "document_ref": "MOCK-0", "title": "mock",
                                   "doc_file_name": "", "doc_mime_type": "", "file_mime_type": "application/json"}],
        "quantities": quantities,
        "data_files": data_files,
        "releases": [{"tag": tag, "release_date": f"2026-0{i + 1}-01T00:00:00+00:00", "comment": "synthetic",
                      "data_files": release_files[tag]} for i, tag in enumerate(MOCK_RELEASES)],
    }
    os.makedirs(MOCK_IMO_DIR, exist_ok=True)
    with open(os.path.join(MOCK_IMO_DIR, "schema.json"), "w") as out:
        json.dump(schema, out, indent=1)


def store_imo_reference(f, workdir):
    imo = lbs.Imo(flatfile_location=MOCK_IMO_DIR)
    rel = "/releases/vMOCK2/satellite"
    g = f.create_group("imo")
    inst = lbs.InstrumentInfo.from_imo(imo, f"{rel}/MOCKFT/instrument_info")
    for key in ("boresight_rotangle_rad", "spin_boresight_angle_rad", "spin_rotangle_rad", "hwp_rpm", "wafer_space_cm"):
        g.attrs[f"instrument_{key}"] = getattr(inst, key)
    write(g, "instrument_bore2spin_quat", inst.bore2spin_quat.quats[0])
    sstr = lbs.SpinningScanningStrategy.from_imo(imo, f"{rel}/scanning_parameters")
    for key in ("spin_sun_angle_rad", "precession_rate_hz", "spin_rate_hz"):
        g.attrs[f"scan_{key}"] = getattr(sstr, key)
    for channel in MOCK_CHANNELS:
        ch = lbs.FreqChannelInfo.from_imo(imo, f"{rel}/MOCKFT/{channel}/channel_info")
        c = g.create_group(f"channel/{channel}")
        for key in ("number_of_detectors", "net_detector_ukrts", "net_channel_ukrts", "bandcenter_ghz", "fmin_hz", "alpha"):
            c.attrs[key] = getattr(ch, key)
        c.attrs["detector_objs"] = [str(x) for x in ch.detector_objs]
        for det_name in MOCK_CHANNELS[channel]:
            det = lbs.DetectorInfo.from_imo(imo, f"{rel}/MOCKFT/{channel}/{det_name}/detector_info")
            d = c.create_group(det_name)
            for key in ("sampling_rate_hz", "fwhm_arcmin", "bandcenter_ghz", "net_ukrts", "fknee_mhz", "pol_angle_rad",
                        "psi_rad", "squid"):
                d.attrs[key] = getattr(det, key)
            for key in ("wafer", "pixel", "pol", "orient", "channel"):
                d.attrs[key] = str(getattr(det, key))
            write(d, "quats", det.quat.quats)
            if det.quat.start_time is not None:
                d.attrs["quat_start_time_s"] = det.quat.start_time
                d.attrs["quat_sampling_rate_hz"] = det.quat.sampling_rate_hz

    # J: observation of the M-100 detectors described entirely by the IMO
    dets = lbs.detector_list_from_parameters(imo, [{"channel_info_obj": f"{rel}/MOCKFT/M-100/channel_info"}])
    sim = lbs.Simulation(base_path=os.path.join(workdir, "J"), start_time=100.0, duration_s=60.0,
                         random_seed=12345, imo=imo)
    sim.set_instrument(inst)
    sim.create_observations(detectors=dets, allocate_tod=False)
    sim.set_scanning_strategy(scanning_strategy=sstr, delta_time_s=1.0 / 19.0)
    sim.set_hwp(lbs.IdealHWP(ang_speed_radpsec=inst.hwp_rpm * 2 * np.pi / 60))
    sim.prepare_pointings()
    obs = sim.observations[0]
    j = f.create_group("J")
    ptg, hwp = obs.get_pointings("all")
    write(j, "pointings_ecl", ptg)
    write(j, "hwp_angle", hwp)
    j.attrs["detector_names"] = [d.name for d in dets]


def main():
    commit = subprocess.run(
        ["git", "-C", os.path.dirname(lbs.__file__), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True,
    ).stdout.strip()
    build_mock_imo()
    with tempfile.TemporaryDirectory() as workdir, h5py.File(OUTPUT, "w") as f:
        f.attrs["litebird_sim_commit"] = commit
        f.attrs["litebird_sim_version"] = lbs.__version__
        f.attrs["numpy_version"] = np.__version__
        f.attrs["astropy_version"] = astropy.__version__
        f.attrs["ducc0_version"] = ducc0.__version__
        f.attrs["stride"] = STRIDE

        c = f.create_group("config")
        for key, value in dict(
            spin_sun_angle_rad=SPIN_SUN_ANGLE_RAD,
            precession_rate_hz=PRECESSION_RATE_HZ,
            spin_rate_hz=SPIN_RATE_HZ,
            spin_boresight_angle_rad=SPIN_BORESIGHT_ANGLE_RAD,
            sampling_rate_hz=SAMPLING_RATE_HZ,
            duration_s=DURATION_S,
            start_time_b=START_TIME_B,
            hwp_ang_speed_radpsec=HWP_ANG_SPEED_RADPSEC,
            hwp_start_angle_rad=HWP_START_ANGLE_RAD,
            spin_rotangle_rad_b=np.pi,
            boresight_rotangle_rad_b=0.3,
        ).items():
            c.attrs[key] = value
        write(c, "detector_angles_deg", np.array(DETECTOR_ANGLES_DEG))
        write(c, "pol_angles_rad", np.array(POL_ANGLES_RAD))
        write(c, "detector_quats_b", np.array([d.quat.quats[0] for d in detectors_b()]))
        write(c, "ecl_to_gal_rot_matrix", ECL_TO_GAL_ROT_MATRIX)
        write(c, "ecl_to_gal_quat", ECL_TO_GAL_QUAT)

        exact = 1.0 / SAMPLING_RATE_HZ

        # A: boresight, no offsets, start at 0, no HWP
        det_a = [lbs.DetectorInfo(name="boresight", sampling_rate_hz=SAMPLING_RATE_HZ)]
        sim = make_sim(workdir, "A", 0.0, instrument_a(), det_a, exact)
        store_pointings(f, "A", sim.observations[0], with_hwp=False)

        # B: three detectors, rotated instrument, absolute start time, HWP
        sim = make_sim(workdir, "B", START_TIME_B, instrument_b(), detectors_b(), exact, hwp=hwp_b())
        store_pointings(f, "B", sim.observations[0], with_hwp=True)
        # E: h-maps of B in both coordinate systems; F: hit map of B
        store_h_maps(f, "E_gal", sim, lbs.CoordinateSystem.Galactic)
        store_h_maps(f, "E_ecl", sim, lbs.CoordinateSystem.Ecliptic)
        g = f.create_group("F")
        write(g, "hitmap", sim.get_hitmap(nside=NSIDE, output_coordinate_system=lbs.CoordinateSystem.Galactic))

        # C: B split into two observations; store the second one
        sim = make_sim(workdir, "C", START_TIME_B, instrument_b(), detectors_b(), exact, hwp=hwp_b(), num_of_obs=2)
        g = store_pointings(f, "C", sim.observations[1], with_hwp=True)
        g.attrs["start_time_simulation_s"] = START_TIME_B

        # D: B with the default anchor spacing (interpolated pointings)
        sim = make_sim(workdir, "D", START_TIME_B, instrument_b(), detectors_b(), 2.0, hwp=hwp_b())
        g = store_pointings(f, "D", sim.observations[0], with_hwp=True)
        g.attrs["delta_time_s"] = 2.0

        # G: time-dependent quaternion for det2, 1 Hz, starting 3.3 s before the observation
        g_start, g_rate, g_n = START_TIME_B - 3.3, 1.0, 605
        nominal = det_quat(*np.deg2rad(DETECTOR_ANGLES_DEG[2]))
        g_quats = wobble_series(g_start, g_n, g_rate, np.deg2rad(1.0 / 60.0), 97.0, nominal)
        for name, delta in (("G_exact", exact), ("G_interp", 2.0)):
            dets = detectors_b()
            dets[2] = lbs.DetectorInfo(
                name="det2", sampling_rate_hz=SAMPLING_RATE_HZ, pol_angle_rad=POL_ANGLES_RAD[2],
                quat=lbs.RotQuaternion(quats=g_quats.copy(), start_time=g_start, sampling_rate_hz=g_rate),
            )
            sim = make_sim(workdir, name, START_TIME_B, instrument_b(), dets, delta, hwp=hwp_b())
            g = store_pointings(f, name, sim.observations[0], with_hwp=False)
            g.attrs["delta_time_s"] = delta
        g = f["G_exact"]
        write(g, "det2_quats", g_quats)
        g.attrs["det2_start_time_s"] = g_start
        g.attrs["det2_sampling_rate_hz"] = g_rate

        # H: spin-to-ecliptic attitude replaced by a 1 Hz time series
        h_series = {}

        def replace_attitude(sim):
            analytic = sim.spin2ecliptic_quats
            quats = np.array(analytic.quats)
            t = START_TIME_B + np.arange(quats.shape[0])
            for k in range(quats.shape[0]):
                lbs.quat_right_multiply(quats[k], *lbs.quat_rotation_y(np.deg2rad(2.0 / 60.0) * np.cos(2 * np.pi * (t[k] - START_TIME_B) / 53.0)))
            h_series["quats"] = quats.copy()
            sim.spin2ecliptic_quats = lbs.RotQuaternion(quats=quats, start_time=START_TIME_B, sampling_rate_hz=1.0)

        sim = make_sim(workdir, "H", START_TIME_B, instrument_b(), detectors_b(), 1.0, hwp=hwp_b(),
                       setup_attitude=replace_attitude)
        g = store_pointings(f, "H", sim.observations[0], with_hwp=False)
        write(g, "attitude_quats", h_series["quats"])
        g.attrs["attitude_start_time_s"] = START_TIME_B
        g.attrs["attitude_sampling_rate_hz"] = 1.0

        # I: analytic attitude times a spin-frame perturbation on the sample grid
        i_series = {}

        def perturb_attitude(sim):
            base = sim.spin2ecliptic_quats
            n = base.quats.shape[0]
            t = START_TIME_B + np.arange(n) * exact
            pert = np.empty((n, 4))
            for k in range(n):
                q = np.array(lbs.quat_rotation_x(np.deg2rad(1.5 / 60.0) * np.sin(2 * np.pi * (t[k] - START_TIME_B) / 41.0)))
                lbs.quat_left_multiply(q, *lbs.quat_rotation_y(np.deg2rad(0.7 / 60.0) * np.cos(2 * np.pi * (t[k] - START_TIME_B) / 23.0)))
                pert[k] = q
            i_series["quats"] = pert.copy()
            sim.spin2ecliptic_quats = base * lbs.RotQuaternion(quats=pert, start_time=START_TIME_B, sampling_rate_hz=SAMPLING_RATE_HZ)

        sim = make_sim(workdir, "I", START_TIME_B, instrument_b(), detectors_b(), exact, hwp=hwp_b(),
                       setup_attitude=perturb_attitude)
        g = store_pointings(f, "I", sim.observations[0], with_hwp=False)
        write(g, "perturbation_quats", i_series["quats"])
        g.attrs["perturbation_start_time_s"] = START_TIME_B
        g.attrs["perturbation_sampling_rate_hz"] = SAMPLING_RATE_HZ

        store_imo_reference(f, workdir)

    print("wrote", OUTPUT, f"({os.path.getsize(OUTPUT) / 2**20:.2f} MiB)")


if __name__ == "__main__":
    main()
