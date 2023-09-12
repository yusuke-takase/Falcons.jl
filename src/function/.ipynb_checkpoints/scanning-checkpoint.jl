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
    @views @inbounds for i = eachindex(ss.quat)
        q_point_idet     = q_point[i]
        q_pol_angle_idet = q_pol_angle[i]
        @views @inbounds for j = eachindex(time_array)
            t              = time_array[j]
            q_revol        = rotate_quat(omega_revol, t, ez)
            q_prec         = rotate_quat(omega_prec,  t, ex)
            q_spin         = rotate_quat(omega_spin,  t, spin_axis)
            
            Q              = q_revol * q_prec * q_spin
            q_point_t      = Q * q_point_idet / Q
            q_pol_t        = Q * q_pol_angle_idet / Q
            
            vec_point      = vect(q_point_t)
            poldir         = vect(q_pol_t)
            
            #θ, ϕ           = vec2ang_minuspi_to_pi(vec_point[1], vec_point[2], vec_point[3])
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


mutable struct ScanningStrategy_ThetaPhi{T<:AbstractFloat, I<:Int, AA<:AbstractArray{T}, AS<:AbstractString}
    nside::I
    duration::I
    sampling_rate::T
    alpha::T
    beta::T
    prec_rpm::T
    spin_rpm::T
    hwp_rpm::T
    FP_theta::AA
    FP_phi::AA
    start_point::AS
    start_angle::T
    coord::AS
end

"""
    gen_ScanningStrategy_ThetaPhi(args***)

This function generate scanning strategy.
    # Arguments
    ...
    nside::Int
    duration::Int
    sampling_rate::Int
    alpha::Float
    beta::Float
    prec_rpm::Float
    spin_rpm::Float
    hwp_rpm::Float
    FP_theta::Array
    FP_phi::Array
    start_point::String
    start_angle::Float
    ...
    
    # Returns 
    ...
    - `scanning_strategy_structure`::ScanningStrategy_ThetaPhi
    ...
"""
function gen_ScanningStrategy_ThetaPhi(;
        nside=128, 
        duration=60*60*24*365, 
        sampling_rate=1.0,
        alpha=45, 
        beta=50, 
        prec_rpm=period2rpm(192.348), 
        spin_rpm=0.05, 
        hwp_rpm=0, 
        FP_theta=[0.0], 
        FP_phi=[0.0], 
        start_point="equator", 
        start_angle=0.0,
        coord="E",
    )
    @warn "This function (gen_ScanningStrategy_ThetaPhi) will not be maintenanced!"
    ScanningStrategy_ThetaPhi(
        nside,
        duration,
        sampling_rate,
        Float64(alpha),
        Float64(beta),
        Float64(prec_rpm),
        Float64(spin_rpm),
        Float64(hwp_rpm),
        Float64.(FP_theta),
        Float64.(FP_phi),
        start_point,
        start_angle,
        coord,
    )
end


"""
    get_pointings_tuple(ss::ScanningStrategy_ThetaPhi, start, stop)

This function will return pointing tod as tuple.
    # Arguments
    ...
    - `ss::ScanningStrategy_ThetaPhi`: the ScanningStrategy_ThetaPhi struct.
    - `start::Int`: the initial time for pointing calculation.
    - `stop::Int`: the finish time for pointing calculation.
    ...

    # Examples
    ```jldoctest
    julia> ss = gen_ScanningStrategy_ThetaPhi()
    julia> pointings = get_pointings_tuple(s, 0, 100)

    # Returns
    ...
    - `pointings[1]`: TOD of theta. shape:(duration*sampling_rate, numOfdet)
    - `pointings[2]`: TOD of phi. shape:(duration*sampling_rate, numOfdet)
    - `pointings[3]`: TOD of psi. shape:(duration*sampling_rate, numOfdet)
    - `pointings[4]`: Array of time. shape:(duration*sampling_rate)
    ...
"""
function get_pointings(ss::ScanningStrategy_ThetaPhi, start, stop)
    @warn "This function (get_pointings(ss::ScanningStrategy_ThetaPhi, start, stop)) will not be maintenanced!"
    resol = Resolution(ss.nside)
    omega_spin = rpm2angfreq(ss.spin_rpm)
    omega_prec = rpm2angfreq(ss.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    time_array = start:1/ss.sampling_rate:stop-1/ss.sampling_rate |> LinRange
    if start > stop-1/ss.sampling_rate
        error("ERROR: \n The `start` time of the calculation is greater than or equal to the `stop` time.")
    end
    loop_times = length(time_array)
    numof_det = length(ss.FP_theta)

    psi_tod = zeros(loop_times, numof_det)
    theta_tod = zeros(loop_times, numof_det)
    phi_tod = zeros(loop_times, numof_det)
    
    ex = @SVector [1.0, 0.0, 0.0]
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    spin_axis = @SVector [cosd(ss.alpha), 0, sind(ss.alpha)]
    qu₀       = Quaternion(0.,0.,-1.,0.)
    qb₀, qd₀ = initial_state(ss)
    
    @views @inbounds for j = eachindex(ss.FP_theta)
        qtheta_in_FP = rotate_quat(deg2rad(ss.FP_theta[j]), ey)
        qphi_in_FP = rotate_quat(deg2rad(ss.FP_phi[j]), vect(qb₀))
        qfp = qphi_in_FP * qtheta_in_FP
        qp₀ⱼ = qfp * qb₀ / qfp
        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            qᵣ = rotate_quat(omega_revol, t, ez)
            qₚ = rotate_quat(omega_prec, t, ex)
            qₛ = rotate_quat(omega_spin, t, spin_axis)
            Q = qᵣ * qₚ * qₛ
            qp = Q * qp₀ⱼ / Q
            qu = Q * qu₀ / Q
            
            p = vect(qp)
            u = vect(qu)
            
            ell = (p × ez) × p  
            θ, ϕ = vec2ang(p[1], p[2], p[3])
            theta_tod[i, j] = θ
            phi_tod[i, j] = ϕ

            k = ell × u
            cosk = dot(u, ell) / (norm(u) * norm(ell))
            cosk = ifelse(abs(cosk) > 1.0, sign(cosk), cosk)
            
            sign_kz = sign(k[3])
            sign_kz = ifelse(sign_kz==0, -1, sign_kz)
            sign_pz = sign(p[3])
            sign_pz = ifelse(sign_pz==0, 1, sign_pz)
            psi_tod[i, j] = acos(cosk) * sign_kz * sign_pz
        end
    end
    if ss.coord == "G"
        rotate_coordinates_e2g!(theta_tod, phi_tod, psi_tod)
    end
    return (theta_tod, phi_tod, psi_tod, time_array)
end


function show_ss(ss::ScanningStrategy_ThetaPhi)
    @printf("%-24s : %i \n", "nside", ss.nside)
    @printf("%-24s : %0.1f \n", "duration [sec]", ss.duration)
    @printf("%-24s : %0.1f \n", "sampling rate [Hz]", ss.sampling_rate)
    @printf("%-24s : %0.1f \n", "alpha [deg]", ss.alpha)
    @printf("%-24s : %0.1f \n", "beta [deg]", ss.beta)
    
    @printf("%-24s : %0.3f\n", "prec. period [min]", period2rpm(ss.prec_rpm))
    @printf("%1s %-22s : %f\n", '\u21B3', "prec. rate [rpm]", ss.prec_rpm)
    
    @printf("%-24s : %0.3f\n", "spin period [min]", period2rpm(ss.spin_rpm))
    @printf("%1s %-22s : %f\n", '\u21B3', "spin rate [rpm]", ss.spin_rpm)
    @printf("%-24s : %f \n", "HWP rot. rate[rpm]", ss.hwp_rpm)
    @printf("%-24s : %s \n", "start point", ss.start_point)
    @printf("%-24s : %f \n", "start angle", ss.start_angle)
    @printf("%-24s : %s \n", "coordinate system", ss.coord)
    @printf("%-24s\n", "FPU")
    for i in eachindex(ss.FP_theta)
        @printf("\u21B3 Det.%i(θ,φ)%-12s : (%0.3f, %0.3f) \n", i, "", ss.FP_theta[i], ss.FP_phi[i])    
    end
end


function initial_state(ss::ScanningStrategy_ThetaPhi)
    @warn "This function (initial_state) will not be maintenanced!"
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    spin_axis = [cosd(ss.alpha), 0, sind(ss.alpha)]
    b₀ = [cosd(ss.alpha+ss.beta), 0, sind(ss.alpha+ss.beta)]
    
    flip_angle = 0
    scan_direction = 1
    if ss.start_point == "equator"
        flip_angle = π
        scan_direction = -1
    end
    
    b₀ = rotate_vec(b₀, flip_angle, spin_axis)
    d₀ = rotate_vec(b₀, -π/2, ey)
    qb₀ = Quaternion(0.0, b₀)
    qd₀ = Quaternion(0.0, d₀)
    q = rotate_quat(ss.start_angle, 1.0, ez)
    qb₀ = q * qb₀ / q
    qd₀ = q * qd₀ / q
    return (qb₀, qd₀)
end