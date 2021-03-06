mutable struct ScanningStrategy{T<:AbstractFloat, I<:Int, AA<:AbstractArray{T}, AS<:AbstractString}
    nside::I
    duration::I
    sampling_rate::I
    alpha::T
    beta::T
    prec_rpm::T
    spin_rpm::T
    hwp_rpm::T
    FP_theta::AA
    FP_phi::AA
    start_point::AS
    start_angle::T
end

rpm2angfreq(rpm) = (2.0π / 60.0) * rpm
#period2rpm(period) = 1.0 / period
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
        nside=128, duration=60*60*24*365, sampling_rate=1, alpha=45, beta=50, 
        prec_rpm=period2rpm(192.348), spin_rpm=0.05, hwp_rpm=0, 
        FP_theta=[0.0], FP_phi=[0.0], start_point="equator", start_angle=0.0)
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
        start_angle
    )
    return scanning_strategy_structure
end

@inline function initial_pointings_rot(ss::ScanningStrategy)
    if ss.start_point =="pole"
        #= Set the initial position of the boresight to near the zenith. =#
        initial_vec = @views @SVector [sind(ss.alpha-ss.beta), 0, cosd(ss.alpha-ss.beta)]
    end
    if ss.start_point == "equator"
        #= Set the initial position of the boresight to near the equator. =#
        initial_vec = @views @SVector [cosd(ss.alpha-ss.beta), 0 , sind(ss.alpha-ss.beta)]
    end
    z_axis = [0, 0, 1]
    bore_quat = Quaternion(0.0, initial_vec)
    
    q = quaternion_rotator(ss.start_angle, 1.0, z_axis)
    bore_quat_rot = q * bore_quat / q
    
    return @views @SVector [bore_quat_rot.q1, bore_quat_rot.q2, bore_quat_rot.q3]
end

@inline function initial_pointings(ss::ScanningStrategy)
    if ss.start_point == "pole"
        #= Set the initial position of the boresight to near the zenith. =#
        initial_vec = @views @SVector [sind(ss.alpha-ss.beta), 0, cosd(ss.alpha-ss.beta)]
    end
    if ss.start_point == "equator"
        #= Set the initial position of the boresight to near the equator. =#
        initial_vec = @views @SVector [cosd(ss.alpha-ss.beta), 0 , sind(ss.alpha-ss.beta)]
    end
    return initial_vec
end

"""
    get_pointings_tuple(SS::ScanningStrategy, start, stop)

This function will return pointing tod as tuple.
    # Arguments
    ...
    - `SS::ScanningStrategy`: the ScanningStrategy struct.
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
@inline function get_pointings_tuple(SS::ScanningStrategy, start, stop)
    resol = Resolution(SS.nside)
    omega_spin = rpm2angfreq(SS.spin_rpm)
    omega_prec = rpm2angfreq(SS.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)

    #time_array = @views Vector(start:1/SS.sampling_rate:stop-1/SS.sampling_rate)
    time_array = start:1/SS.sampling_rate:stop-1/SS.sampling_rate
    loop_times = length(time_array)
    numof_det = length(SS.FP_theta)

    psi_tod = zeros(loop_times, numof_det)
    theta_tod = zeros(loop_times, numof_det)
    phi_tod = zeros(loop_times, numof_det)
    
    antisun_axis = @SVector [1.0, 0.0, 0.0]
    y_axis = @SVector [0.0, 1.0, 0.0]
    z_axis = @SVector [0.0, 0.0, 1.0]
    
    antisun_axis = vector_rotator(antisun_axis, SS.start_angle, z_axis)
    spin_axis = @SVector [cosd(SS.alpha), 0.0, sind(SS.alpha)]
    spin_axis = vector_rotator(spin_axis, SS.start_angle, z_axis)
    
    #bore_0 = initial_pointings(SS)
    bore_0 = initial_pointings_rot(SS)

    boresight_0 = Quaternion(0.0, bore_0)
    #detector_vec_0 = @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_vec_0 = @SVector [-boresight_0.q3, boresight_0.q2, boresight_0.q1]
    detector_orientation_0 = Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = bore_0 × detector_vec_0
    travel_direction_0 = Quaternion(0.0, travel_direction_vec_0)

    @views @inbounds @simd for j = eachindex(SS.FP_theta)
        q_theta_in_FP = quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis)
        q_phi_in_FP = quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP

        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            q_revol = quaternion_rotator(omega_revol, t, z_axis)
            q_prec = quaternion_rotator(omega_prec, t, antisun_axis)
            q_spin = quaternion_rotator(omega_spin, t, spin_axis)
            Q = q_revol * q_prec * q_spin
            P = Q * pointing / Q
            q_det = Q * detector_orientation_0 / Q
            travel_direction = Q * travel_direction_0 / Q
            pointing_t = @SVector [P.q1, P.q2, P.q3]
            travel_direction_vec = @SVector [travel_direction.q1, travel_direction.q2, travel_direction.q3]

            longitude = pointing_t × (pointing_t × z_axis)
            ang_t = vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            theta_tod[i, j] = ang_t[1]
            phi_tod[i, j] = ang_t[2]

            divergent_vec = longitude × travel_direction_vec
            cosK = dot(travel_direction_vec, longitude) / (norm(travel_direction_vec) * norm(longitude))

            cosK = ifelse(abs(cosK) > 1.0, sign(cosK), cosK)

            divergent_vec3 = sign(divergent_vec[3])
            divergent_vec3 = ifelse(divergent_vec3==0, -1, divergent_vec3)
            pointing_t3 = sign(pointing_t[3])
            pointing_t3 = ifelse(pointing_t3==0, 1, pointing_t3)

            psi_tod[i, j] = acos(cosK) * divergent_vec3 * pointing_t3
        end
    end
    return (theta_tod, phi_tod, psi_tod, time_array)
end

@inline function get_pointings_xyz_tuple(SS::ScanningStrategy, start, stop)
    resol = Resolution(SS.nside)
    omega_spin = rpm2angfreq(SS.spin_rpm)
    omega_prec = rpm2angfreq(SS.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)

    #time_array = @views Vector(start:1/SS.sampling_rate:stop-1/SS.sampling_rate)
    time_array = start:1/SS.sampling_rate:stop-1/SS.sampling_rate
    loop_times = length(time_array)
    numof_det = length(SS.FP_theta)

    psi_tod = zeros(loop_times, numof_det)
    #theta_tod = zeros(loop_times, numof_det)
    #phi_tod = zeros(loop_times, numof_det)
    xyz_tod = zeros(3, loop_times, numof_det)
    
    antisun_axis = @SVector [1.0, 0.0, 0.0]
    y_axis = @SVector [0.0, 1.0, 0.0]
    z_axis = @SVector [0.0, 0.0, 1.0]
    
    antisun_axis = vector_rotator(antisun_axis, SS.start_angle, z_axis)
    spin_axis = @SVector [cosd(SS.alpha), 0.0, sind(SS.alpha)]
    spin_axis = vector_rotator(spin_axis, SS.start_angle, z_axis)
    
    #bore_0 = initial_pointings(SS)
    bore_0 = initial_pointings_rot(SS)

    boresight_0 = Quaternion(0.0, bore_0)
    #detector_vec_0 = @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_vec_0 = @SVector [-boresight_0.q3, boresight_0.q2, boresight_0.q1]
    detector_orientation_0 = Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = bore_0 × detector_vec_0
    travel_direction_0 = Quaternion(0.0, travel_direction_vec_0)

    @views @inbounds @simd for j = eachindex(SS.FP_theta)
        q_theta_in_FP = quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis)
        q_phi_in_FP = quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP

        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            q_revol = quaternion_rotator(omega_revol, t, z_axis)
            q_prec = quaternion_rotator(omega_prec, t, antisun_axis)
            q_spin = quaternion_rotator(omega_spin, t, spin_axis)
            Q = q_revol * q_prec * q_spin
            P = Q * pointing / Q
            q_det = Q * detector_orientation_0 / Q
            travel_direction = Q * travel_direction_0 / Q
            pointing_t = @SVector [P.q1, P.q2, P.q3]
            travel_direction_vec = @SVector [travel_direction.q1, travel_direction.q2, travel_direction.q3]

            longitude = pointing_t × (pointing_t × z_axis)
            ang_t = vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            #theta_tod[i, j] = ang_t[1]
            #phi_tod[i, j] = ang_t[2]
            xyz_tod[:, i, j] = pointing_t

            divergent_vec = longitude × travel_direction_vec
            cosK = dot(travel_direction_vec, longitude) / (norm(travel_direction_vec) * norm(longitude))

            cosK = ifelse(abs(cosK) > 1.0, sign(cosK), cosK)

            divergent_vec3 = sign(divergent_vec[3])
            divergent_vec3 = ifelse(divergent_vec3==0, -1, divergent_vec3)
            pointing_t3 = sign(pointing_t[3])
            pointing_t3 = ifelse(pointing_t3==0, 1, pointing_t3)

            psi_tod[i, j] = acos(cosK) * divergent_vec3 * pointing_t3
        end
    end
    return (xyz_tod, psi_tod, time_array)
end

"""
    get_pointing_pixels(SS::ScanningStrategy, start, stop)

This function will return pointing pixel tod as tuple.
    # Arguments
    ...
    - `SS::ScanningStrategy`: the ScanningStrategy struct.
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
@inline function get_pointing_pixels(SS::ScanningStrategy, start, stop)
    theta_tod, phi_tod, psi_tod, time_array = get_pointings_tuple(SS, start, stop)
    loop_times = length(theta_tod[:,1])
    numOfdet = length(theta_tod[1,:])
    resol = Resolution(SS.nside)

    pix_tod = @views zeros(Int64, loop_times, numOfdet)
    for j in eachindex(theta_tod[1,:])
        for i in eachindex(theta_tod[:,1])
            pix_tod[i, j] = @views ang2pixRing(resol, theta_tod[i, j], phi_tod[i, j])
        end
    end
    return (pix_tod, psi_tod, time_array)
end

"""
    get_pointings(SS::ScanningStrategy, start, stop)

This function will return pointing tod as dictionary type.
    # Arguments
    ...
    - `SS::ScanningStrategy`: the ScanningStrategy struct.
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
function get_pointings(SS::ScanningStrategy, start, stop)
    pointings = get_pointings_tuple(SS, start, stop)
    pointings = Dict{String,AbstractArray{Float64}}(
        "theta" => pointings[1],
        "phi" => pointings[2],
        "psi" => pointings[3],
        "time" => pointings[4]
    )
    return pointings
end
