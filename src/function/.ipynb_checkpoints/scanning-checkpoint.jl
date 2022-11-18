@inline function vec2ang_ver2(x, y, z)
    norm = sqrt(x^2 + y^2 + z^2)
    theta = acos(z / norm)
    phi = atan(y, x)
    phi = ifelse(phi > π, π-phi, phi)
    return (theta, phi)
end


@inline function quaternion_rotator(omega, t, rotate_axis)
    #= Generate a quaternion that rotates by the angle omega*t around the rotate_axis axis. =#
    rot_ang = omega * t
    Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
end

@inline function vector_rotator(vec, rot_ang, rotate_axis)
    q = Quaternion([cos(rot_ang/2.0), rotate_axis[1]*sin(rot_ang/2.0), rotate_axis[2]*sin(rot_ang/2.0), rotate_axis[3]*sin(rot_ang/2.0)])
    vec_q = Quaternion(0.0, vec)
    rot_vec = q * vec_q / q
    return vect(rot_vec)
end

# θ, φ in ecliptic coord → θ, φ in galactic coord 
function rot_E2G_ang(β, λ)
    # reference:: https://aas.aanda.org/articles/aas/full/1998/01/ds1449/node3.html
    # equinox 1950
    # (theta,phi)
    # -π < λ,l < π  ,   0 < β < π
    β_NGP = deg2rad(29.81)
    λ_0 = deg2rad(269.32)
    l_1 = deg2rad(6.38)
    β = pi/2.0 - β # 0 < β < π →  -π/2 < β < π/2
    sin_b = (sin(β) * sin(β_NGP) - cos(β) * cos(β_NGP) * sin(λ - λ_0))
    b = asin(sin_b)

    cos_l_l_1 = (cos(λ - λ_0) * cos(β) / cos(b))
    sin_l_l_1 = ((sin(β) * cos(β_NGP) + cos(β) * sin(β_NGP) * sin(λ - λ_0)) / cos(b))
    l_l_1 = atan(sin_l_l_1, cos_l_l_1)
    l = l_l_1 + l_1
    b = pi/2.0 - b
    return b, l
end

# vector in ecliptic coord → vector in galactic coord. But output vector is static vector
function ecliptic2galactic(v)
    β, λ = vec2ang(v[1], v[2], v[3])
    b, l = rot_E2G_ang(β, λ)
    x, y, z = ang2vec(b, l)
    return @SVector [x, y, z]
end


mutable struct ScanningStrategy{T<:AbstractFloat, I<:Int, AA<:AbstractArray{T}, AS<:AbstractString}
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

rpm2angfreq(rpm) = (2.0π / 60.0) * rpm

function period2rpm(period,; unit="min")
    if unit == "min"
        rpm = 1.0 / period
    end
    if unit == "sec"
        rpm = 1.0 / (period/60.0)
    end
    if unit == "hour"
        rpm = 1.0 / (period/60.0/60.0)
    end
    return rpm
end

function rpm2period(rpm,; unit="min")
    if unit == "min"
        period = 1.0 / rpm
    end
    if unit == "sec"
        period = 1.0 / (rpm/60.0)
    end
    if unit == "hour"
        period = 1.0 / (rpm/60.0/60.0)
    end
    return period
end
    
function show_ss(ss::ScanningStrategy)
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
    @printf("%-24s\n", "FPU")
    for i in eachindex(ss.FP_theta)
        @printf("\u21B3 Det.%i(θ,φ)%-12s : (%0.3f, %0.3f) \n", i, "", ss.FP_theta[i], ss.FP_phi[i])    
    end
end

"""
    gen_ScanningStrategy(args***)

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
    - `scanning_strategy_structure`::ScanningStrategy
    ...
"""
function gen_ScanningStrategy(;
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
        coord="E" 
    )
    scanning_strategy_structure = ScanningStrategy(
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
        coord
    )
    return scanning_strategy_structure
end

function initial_state(ss::ScanningStrategy)
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    spin_axis = [cosd(ss.alpha), 0, sind(ss.alpha)]
    b₀ = [cosd(ss.alpha+ss.beta), 0, sind(ss.alpha+ss.beta)]
    anti_sun_axis = ex
    
    flip_angle = 0
    scan_direction = 1
    if ss.start_point == "equator"
        flip_angle = π
        scan_direction = -1
    end
    
    b₀ = vector_rotator(b₀, flip_angle, spin_axis)
    d₀ = vector_rotator(b₀, -π/2, ey)
    u₀ = scan_direction *  b₀ × d₀
    
    if ss.coord == "G"
        ex = ecliptic2galactic(ex)
        ey = ecliptic2galactic(ey)
        ez = ecliptic2galactic(ez)
        spin_axis = ecliptic2galactic(spin_axis)
        b₀ = ecliptic2galactic(b₀)
        anti_sun_axis = ex
        d₀ = ecliptic2galactic(d₀)
        u₀ = ecliptic2galactic(u₀)
    end
    
    qb₀ = Quaternion(0.0, b₀)
    qd₀ = Quaternion(0.0, d₀)
    qu₀ = Quaternion(0.0, u₀)
    qspin_axis = Quaternion(0.0, spin_axis)
    qanti_sun_axis = Quaternion(0.0, anti_sun_axis)

    q = quaternion_rotator(ss.start_angle, 1.0, ez)
    qb₀ = q * qb₀ / q
    qd₀ = q * qd₀ / q
    qu₀ = q * qu₀ / q
    qspin_axis = q * qspin_axis / q
    qanti_sun_axis = q * qanti_sun_axis / q

    spin_axis = vect(qspin_axis)
    anti_sun_axis = vect(qanti_sun_axis)
    
    return (qb₀, qd₀, qu₀, spin_axis, anti_sun_axis)
end

"""
    get_pointings_tuple(ss::ScanningStrategy, start, stop)

This function will return pointing tod as tuple.
    # Arguments
    ...
    - `ss::ScanningStrategy`: the ScanningStrategy struct.
    - `start::Int`: the initial time for pointing calculation.
    - `stop::Int`: the finish time for pointing calculation.
    ...

    # Examples
    ```jldoctest
    julia> ss = gen_ScanningStrategy()
    julia> pointings = get_pointings_tuple(s, 0, 100)

    # Returns
    ...
    - `pointings[1]`: TOD of theta. shape:(duration*sampling_rate, numOfdet)
    - `pointings[2]`: TOD of phi. shape:(duration*sampling_rate, numOfdet)
    - `pointings[3]`: TOD of psi. shape:(duration*sampling_rate, numOfdet)
    - `pointings[4]`: Array of time. shape:(duration*sampling_rate)
    ...
"""
function get_pointings_tuple(ss::ScanningStrategy, start, stop)
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
    
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    if ss.coord == "G"
        ey = ecliptic2galactic(ey)
        ez = ecliptic2galactic(ez)
    end

    qb₀, qd₀, qu₀, spin_axis, antisun_axis = initial_state(ss)
    
    @views @inbounds for j = eachindex(ss.FP_theta)
        qtheta_in_FP = quaternion_rotator(deg2rad(ss.FP_theta[j]), 1, ey)
        qphi_in_FP = quaternion_rotator(deg2rad(ss.FP_phi[j]), 1, vect(qb₀))
        qfp = qphi_in_FP * qtheta_in_FP
        qp₀ⱼ = qfp * qb₀ / qfp
        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            qᵣ = quaternion_rotator(omega_revol, t, ez)
            qₚ = quaternion_rotator(omega_prec, t, antisun_axis)
            qₛ = quaternion_rotator(omega_spin, t, spin_axis)
            Q = qᵣ * qₚ * qₛ
            qp = Q * qp₀ⱼ / Q
            qu = Q * qu₀ / Q
            
            p = vect(qp)
            u = vect(qu)
            
            ell = (p × ez) × p  
            
            θ, ϕ = vec2ang_ver2(p[1], p[2], p[3])
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
@inline function get_pointing_pixels(ss::ScanningStrategy, start, stop)
    theta_tod, phi_tod, psi_tod, time_array = get_pointings_tuple(ss, start, stop)
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

"""
    get_pointings(ss::ScanningStrategy, start, stop)

This function will return pointing tod as dictionary type.
    # Arguments
    ...
    - `ss::ScanningStrategy`: the ScanningStrategy struct.
    - `start::Int`: the initial time for pointing calculation.
    - `stop::Int`: the finish time for pointing calculation.
    ...

    # Examples
    ```jldoctest
    julia> ss = gen_ScanningStrategy()
    julia> pointings = get_pointings(s, 0, 100)

    # Returns
    ...
    - `pointings["theta"]`: TOD of theta. shape:(duration*sampling_rate, numOfdet)
    - `pointings["phi"]`: TOD of phi. shape:(duration*sampling_rate, numOfdet)
    - `pointings["psi"]`: TOD of psi. shape:(duration*sampling_rate, numOfdet)
    - `pointings["time"]`: Array of time. shape:(duration*sampling_rate)
    ...
"""
function get_pointings(ss::ScanningStrategy, start, stop)
    pointings = get_pointings_tuple(ss, start, stop)
    pointings = Dict{String,AbstractArray{Float64}}(
        "theta" => pointings[1],
        "phi" => pointings[2],
        "psi" => pointings[3],
        "time" => pointings[4]
    )
    return pointings
end

function ScanningStrategy2map(ss,; division, nmax=5)
    outmap = get_hn_map(ss, division=division, nmax=nmax)
    hmap = outmap["h"] .|> abs2
    return (outmap["hitmap"], hmap)
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


function convert_maps(healpy_maps)
    PolarizedHealpixMap{Float64,RingOrder}(
        healpy_maps[1,:], 
        healpy_maps[2,:], 
        healpy_maps[3,:],
    )
end


