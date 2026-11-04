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
"""

import os
import subprocess
import tempfile

import astropy
import ducc0
import h5py
import numpy as np

import litebird_sim as lbs
from litebird_sim.coordinates import ECL_TO_GAL_QUAT, ECL_TO_GAL_ROT_MATRIX

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "lbs_reference.h5")

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


def main():
    commit = subprocess.run(
        ["git", "-C", os.path.dirname(lbs.__file__), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True,
    ).stdout.strip()
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

    print("wrote", OUTPUT, f"({os.path.getsize(OUTPUT) / 2**20:.2f} MiB)")


if __name__ == "__main__":
    main()
