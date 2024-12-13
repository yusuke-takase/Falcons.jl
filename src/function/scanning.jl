@inline function vec2ang_minuspi_to_pi(x, y, z)
    norm = sqrt(x^2 + y^2 + z^2)
    theta = acos(z / norm)
    phi = atan(y, x)
    phi = ifelse(phi > π, π-phi, phi)
    return (theta, phi)
end

@inline function rotate_quat(omega, t, rotate_axis)
    #= Generate a quaternion that rotates by the angle omega*t around the rotate_axis axis. =#
    rot_ang = mod2pi(omega * t)
    Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
end

@inline function rotate_quat(ϕ, rotate_axis)
    Quaternion([cos(ϕ/2.), rotate_axis[1]*sin(ϕ/2.), rotate_axis[2]*sin(ϕ/2.), rotate_axis[3]*sin(ϕ/2.)])
end

@inline function rotate_vec(vec, rot_ang, rotate_axis)
    q = Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
    vec_q = Quaternion(0.0, vec)
    rot_vec = q * vec_q / q
    return vect(rot_vec)
end
mutable struct ScanningStrategy{T<:AbstractFloat, I<:Int, AS<:AbstractString}
    nside::I
    duration::I
    sampling_rate::T
    alpha::T
    beta::T
    gamma::T
    prec_rpm::T
    spin_rpm::T
    hwp_rpm::T
    start_point::AS
    start_angle::T
    coord::AS
    quat::Vector{Vector{Float64}}
    name::Vector{String}
    info::DataFrame
end

function show_ss(ss::ScanningStrategy)
    @printf("%-24s : %i \n", "nside", ss.nside)
    @printf("%-24s : %0.1f \n", "duration [sec]", ss.duration)
    @printf("%-24s : %0.1f \n", "sampling rate [Hz]", ss.sampling_rate)
    @printf("%-24s : %0.1f \n", "alpha [deg]", ss.alpha)
    @printf("%-24s : %0.1f \n", "beta [deg]", ss.beta)
    @printf("%-24s : %0.1f \n", "gamma [deg]", ss.gamma)

    @printf("%-24s : %0.3f\n", "prec. period [min]", period2rpm(ss.prec_rpm))
    @printf("%1s %-22s : %f\n", '\u21B3', "prec. rate [rpm]", ss.prec_rpm)

    @printf("%-24s : %0.3f\n", "spin period [min]", period2rpm(ss.spin_rpm))
    @printf("%1s %-22s : %f\n", '\u21B3', "spin rate [rpm]", ss.spin_rpm)
    @printf("%-24s : %f \n", "HWP rot. rate[rpm]", ss.hwp_rpm)
    #@printf("%-24s : %s \n", "start point", ss.start_point)
    @printf("%-24s : %f \n", "start angle", ss.start_angle)
    @printf("%-24s : %s \n", "coordinate system", ss.coord)
    @printf("%-24s\n", "FPU")
    for i in eachindex(ss.quat)
        @printf("\u21B3 Det. %i  name | %s %-12s : (x,y,z,w) = [%0.3f, %0.3f, %0.3f, %0.3f] \n", i, ss.name[i], "",
            ss.quat[i][1],
            ss.quat[i][2],
            ss.quat[i][3],
            ss.quat[i][4]
        )
    end
end

mutable struct Pointings
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
    θ::AbstractFloat # position of sky
    φ::AbstractFloat # position of sky
    Ω::Int           # Sky pixel index
    ψ::AbstractFloat # Crossing angle
    ϕ::AbstractFloat # HWP angle
    ξ::AbstractFloat # mod2pi(2ϕ) + ψ
end

function pointings(resol::Resolution, θ, φ, ψ, ϕ)
    vec = ang2vec(θ, φ)
    Ω = ang2pixRing(resol, θ, φ)
    ξ = mod2pi(2ϕ) + ψ
    return Pointings(vec[1],vec[2],vec[3], θ, φ, Ω, ψ, ϕ, ξ)
end

function gen_ScanningStrategy(;
        nside=128,
        duration=60*60*24*365,
        sampling_rate=1.0,
        alpha=45,
        beta=50,
        gamma=0,
        prec_rpm=period2rpm(192.348),
        spin_rpm=0.05,
        hwp_rpm=0,
        start_point="equator",
        start_angle=0.0,
        coord="E",
        quat=[[0.,0.,1.,0.]],
        name= ["boresight"],
        info = DataFrame()
    )
    ScanningStrategy(
        nside,
        duration,
        sampling_rate,
        Float64(alpha),
        Float64(beta),
        Float64(gamma),
        Float64(prec_rpm),
        Float64(spin_rpm),
        Float64(hwp_rpm),
        start_point,
        start_angle,
        coord,
        quat,
        name,
        info
    )
end

function get_satellite(satellite)
    six_month_min = 24*60*60*(365.25/2)/60
    four_day_min  = 24*60*60*4/60
    satellites = ["WMAP", "Planck",      "PICO", "CORE",       "EPIC" , "LiteBIRD"]
    Alpha      = [70.,    7.5,           26,     30,            45    , 45        ]
    Beta       = [22.5,   85,            69,     65,            55    , 50        ]
    T_alpha    = [60.,    six_month_min, 10*60,  four_day_min,  3.2*60, 192.348   ] # min
    T_beta     = [129/60, 1,             1,      2,             1     , 20.       ] # min

    if satellite in satellites
        idx = findfirst(x -> x == satellite, satellites)
    else
        return "Satellite not found in the list."
    end
    ss = gen_ScanningStrategy(;
            nside=128,
            duration=60*60*24*365,
            sampling_rate=1.0,
            alpha=Alpha[idx],
            beta=Beta[idx],
            gamma=0,
            prec_rpm=period2rpm(T_alpha[idx]),
            spin_rpm=period2rpm(T_beta[idx]),
            hwp_rpm=0,
            start_point="equator",
            start_angle=0.0,
            coord="E",
            quat=[[0.,0.,1.,0.]],
            name= ["boresight"],
            info = DataFrame()
    )
end


function imo2ecl_coordinates(ss::ScanningStrategy)
    #=
    This function is loading imo directory. The imo v2.0< have a bug that LFT orientation.
    In order to modify the bug, the comment of `g_gamma` have to be erase.
    =#
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    q_boresight = Quaternion(0.,0.,0.,1.)
    q_pol       = Quaternion(0.,1.,0.,0.)
    q_dets      = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_pol_dets  = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_gamma     = 0.
    polang      = 0.
    for i in eachindex(ss.quat)
        if ss.name[i] == "boresight"
            if ss.start_point == "pole"
                flip = π
            elseif ss.start_point == "equator"
                flip = 0.
            end
            Q              = rotate_quat(deg2rad(90.0-ss.alpha), ey) * rotate_quat(flip, ez) * rotate_quat(deg2rad(ss.beta), ey)
            q_dets[i]      = Q * q_boresight / Q
            q_pol_dets[i]  = Q * q_pol / Q
        else
            telescope  = split.(ss.name[i], "_")[1]
            q_imo      = Quaternion(ss.quat[i][4], ss.quat[i][1], ss.quat[i][2], ss.quat[i][3])
            if telescope     == "000" #LFT
                #q_gamma = rotate_quat(deg2rad(270), ez)
                q_gamma = rotate_quat(deg2rad(ss.gamma), ez)
            elseif telescope == "001" #MFT
                #q_gamma = rotate_quat(deg2rad(240), ez)
                q_gamma = rotate_quat(deg2rad(ss.gamma), ez)
            elseif telescope == "002" #HFT
                #q_gamma = rotate_quat(deg2rad(30), ez)
                q_gamma = rotate_quat(deg2rad(ss.gamma), ez)
            end
            Q = rotate_quat(deg2rad(ss.beta), ey) * q_gamma * q_imo
            if telescope == "001" #MFT
                Q = rotate_quat(π, ez) * Q
            elseif telescope == "002" #HFT
                Q = rotate_quat(π, ez) * Q
            end
            Q              = rotate_quat(deg2rad(90.0-ss.alpha), ey) * Q
            q_dets[i]      = Q * q_boresight / Q
            q_pol_dets[i]  = Q * q_pol / Q
            #q_scan_dets[i] = Q * q_scan / Q
        end
    end
    return (q_dets, q_pol_dets)
end

function _clip_sincos(x)
    # np.clip would be perfect here, but it is slightly inefficient to use here
    # because `x` might either be a float64 or a float32, and we would need to
    # define the −1 and 1 constant accordingly. Better to resort to a combination
    # of calls to min/max, as they are efficiently optimized by LLVM.
    return min(max(x, -1), 1)
end

function polarization_angle(θ, ϕ, poldir)
    cos_psi = _clip_sincos(-sin(ϕ) * poldir[1] + cos(ϕ) * poldir[2])
    sin_psi = _clip_sincos(
          (-cos(θ) * cos(ϕ) * poldir[1])
        + (-cos(θ) * sin(ϕ) * poldir[2])
        + (sin(θ)  * poldir[3])
    )
    return atan(sin_psi, cos_psi)
end

function get_pointings(ss::ScanningStrategy, start, stop)
    resol = Resolution(ss.nside)
    omega_spin = rpm2angfreq(ss.spin_rpm)
    omega_prec = rpm2angfreq(ss.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    time_array = start:1/ss.sampling_rate:stop-1/ss.sampling_rate |> LinRange
    if start > stop-1/ss.sampling_rate
        error("ERROR: \n The `start` time of the calculation is greater than or equal to the `stop` time.")
    end
    loop_times = length(time_array)
    numof_det = length(ss.quat)

    psi_tod = zeros(loop_times, numof_det)
    theta_tod = zeros(loop_times, numof_det)
    phi_tod = zeros(loop_times, numof_det)

    ex = @SVector [1.0, 0.0, 0.0]
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    spin_axis = @SVector [cosd(ss.alpha), 0, sind(ss.alpha)]
    q_scan_direction     = Quaternion(0.,0.,1.,0.)
    q_point, q_pol_angle = imo2ecl_coordinates(ss)
    q_solar_system       = rotate_quat(ss.start_angle, ez)
    @views @inbounds for i = eachindex(ss.quat)
        q_point_idet     = q_point[i]
        q_pol_angle_idet = q_pol_angle[i]
        @views @inbounds for j = eachindex(time_array)
            t              = time_array[j]
            q_revol        = rotate_quat(omega_revol, t, ez)
            q_prec         = rotate_quat(omega_prec,  t, ex)
            q_spin         = rotate_quat(omega_spin,  t, spin_axis)

            Q              = q_solar_system * q_revol * q_prec * q_spin
            q_point_t      = Q * q_point_idet / Q
            q_pol_t        = Q * q_pol_angle_idet / Q

            vec_point      = vect(q_point_t)
            poldir         = vect(q_pol_t)

            θ, ϕ           = vec2ang(vec_point[1], vec_point[2], vec_point[3])
            theta_tod[j,i] = θ
            phi_tod[j,i]   = ϕ
            psi_tod[j,i]   = polarization_angle(θ, ϕ, poldir)
        end
    end
    if ss.coord == "G"
        rotate_coordinates_e2g!(theta_tod, phi_tod, psi_tod)
    end
    return (theta_tod, phi_tod, psi_tod, time_array)
end

function get_pointings(ss::ScanningStrategy, start, stop, step)
    resol = Resolution(ss.nside)
    omega_spin = rpm2angfreq(ss.spin_rpm)
    omega_prec = rpm2angfreq(ss.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    time_array = start:1/ss.sampling_rate:stop-1/ss.sampling_rate |> LinRange
    if start > stop-1/ss.sampling_rate
        error("ERROR: \n The `start` time of the calculation is greater than or equal to the `stop` time.")
    end
    loop_times = length(time_array)
    numof_det = length(ss.quat)

    psi_tod = zeros(loop_times, numof_det)
    theta_tod = zeros(loop_times, numof_det)
    phi_tod = zeros(loop_times, numof_det)

    ex = @SVector [1.0, 0.0, 0.0]
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    spin_axis = @SVector [cosd(ss.alpha), 0, sind(ss.alpha)]
    q_scan_direction     = Quaternion(0.,0.,1.,0.)
    q_point, q_pol_angle = imo2ecl_coordinates(ss)
    q_solar_system       = rotate_quat(ss.start_angle, ez)
    deg_per_day = 360/365.25
    @views @inbounds for i = eachindex(ss.quat)
        q_point_idet     = q_point[i]
        q_pol_angle_idet = q_pol_angle[i]
        @views @inbounds for j = eachindex(time_array)
            t              = time_array[j]
            revol_angle    = deg2rad(deg_per_day * floor(t/step))
            q_revol        = rotate_quat(revol_angle, ez)
            q_prec         = rotate_quat(omega_prec,  t, ex)
            q_spin         = rotate_quat(omega_spin,  t, spin_axis)

            Q              = q_solar_system * q_revol * q_prec * q_spin
            q_point_t      = Q * q_point_idet / Q
            q_pol_t        = Q * q_pol_angle_idet / Q

            vec_point      = vect(q_point_t)
            poldir         = vect(q_pol_t)

            θ, ϕ           = vec2ang(vec_point[1], vec_point[2], vec_point[3])
            theta_tod[j,i] = θ
            phi_tod[j,i]   = ϕ
            psi_tod[j,i]   = polarization_angle(θ, ϕ, poldir)
        end
    end
    if ss.coord == "G"
        rotate_coordinates_e2g!(theta_tod, phi_tod, psi_tod)
    end
    return (theta_tod, phi_tod, psi_tod, time_array)
end

"""
    get_pointing_pixels(ss::ScanningStrategy, start, stop)

This function will return pointing pixel tod as tuple.
    # Arguments
    ...
    - `ss::ScanningStrategy`: the ScanningStrategy struct.
    - `start::Int`: the initial time for pointing calculation.
    - `stop::Int`: the finish time for pointing calculation.
    ...

    # Examples
    ```jldoctest
    julia> ss = gen_ScanningStrategy()
    julia> pointings = get_pointing_pixels(s, 0, 100)

    # Returns
    ...
    - `pointings[1]`: TOD of pixel indicies. shape:(duration*sampling_rate, numOfdet)
    - `pointings[2]`: TOD of psi. shape:(duration*sampling_rate, numOfdet)
    - `pointings[3]`: Array of time. shape:(duration*sampling_rate)
    ...
"""
@inline function get_pointing_pixels(ss, start, stop)
    theta_tod, phi_tod, psi_tod, time_array = get_pointings(ss, start, stop)
    loop_times = length(theta_tod[:,1])
    numOfdet = length(theta_tod[1,:])
    resol = Resolution(ss.nside)

    pix_tod = @views zeros(Int64, loop_times, numOfdet)
    for j in eachindex(theta_tod[1,:])
        for i in eachindex(theta_tod[:,1])
            pix_tod[i, j] = @views ang2pixRing(resol, theta_tod[i, j], phi_tod[i, j])
        end
    end
    return (pix_tod, psi_tod, time_array)
end

function angtod2hitmap(nside::Int, theta_tod, phi_tod)
    resol = Resolution(nside)
    hit_map = zeros(resol.numOfPixels)
    for j in eachindex(theta_tod[1,:])
        theta_tod_jth_det = @views theta_tod[:,j]
        phi_tod_jth_det = @views phi_tod[:,j]
        @inbounds for k in eachindex(theta_tod[:,1])
            ipix = @views ang2pixRing(resol, theta_tod_jth_det[k], phi_tod_jth_det[k])
            hit_map[ipix] += 1
        end
    end
    return hit_map
end

function pixtod2hitmap(nside::Int, pixtod)
    resol = Resolution(nside)
    hit_map = zeros(resol.numOfPixels)
    for i in eachindex(pixtod)
        hit_map[pixtod[i]] += 1
    end
    return hit_map
end

function normarize!(resol::Resolution, maps::Array, hitmap::Array)
    if size(maps) == (3,1,resol.numOfPixels)
        for i in 1:3
            maps[i,1,:] ./= hitmap
        end
    end
    if size(maps) == (3,3,resol.numOfPixels)
        for i in 1:3
            for j in 1:3
                maps[i,j,:] ./= hitmap
            end
        end
    end
    return maps
end
