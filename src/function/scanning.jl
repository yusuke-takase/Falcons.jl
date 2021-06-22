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
end

rpm2angfreq(rpm) = (2.0π / 60.0) * rpm
period2rpm(period) = 1.0 / period

function gen_ScanningStrategy(;nside=128, duration=60*60*24*365, sampling_rate=1, alpha=45, beta=50, prec_rpm=0.005, spin_rpm=0.01, hwp_rpm=0, FP_theta=[0.0], FP_phi=[0.0], start_point="equator")
    scanning_strategy_structure = ScanningStrategy(nside,
        duration,
        sampling_rate,
        Float64(alpha), 
        Float64(beta),
        Float64(prec_rpm),
        Float64(spin_rpm),
        Float64(hwp_rpm),
        Float64.(FP_theta),
        Float64.(FP_phi),
        start_point
    )
    return scanning_strategy_structure
end

@inline function initial_pointings(ss::ScanningStrategy)
    if ss.start_point =="pole"
        #= Set the initial position of the boresight to near the zenith. =#
        initial_vec = @views @SVector [sind(ss.alpha-ss.beta), 0, cosd(ss.alpha-ss.beta)]
    end
    if ss.start_point == "equator"
        #= Set the initial position of the boresight to near the equator. =#
        initial_vec = @views @SVector [cosd(ss.alpha-ss.beta), 0 , sind(ss.alpha-ss.beta)]
    end
    initial_vec
end

#=
@inline function rotate_by_quarternion!(target_quart, quart1, quart2)
    #q_theta_in_FP = quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis)
    #q_phi_in_FP = quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0)
    #q_for_FP = quart1 * quart1
    return quart2 * quart1 * target_quart / quart1 / quart2
end
=#

@inline function get_pointings(SS::ScanningStrategy, start, stop)
    resol = Resolution(SS.nside)

    omega_spin = rpm2angfreq(SS.spin_rpm)
    omega_prec = rpm2angfreq(SS.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    
    #= Compute the TOD of the hitpixel and the TOD of the detector orientation at a specified sampling rate from time start to stop. =#
    time_array = @views Vector(start:1/SS.sampling_rate:stop-1/SS.sampling_rate)
    loop_times = length(time_array)
    numof_det = length(SS.FP_theta)
    
    psi_tod = @views zeros(loop_times, numof_det)
    theta_tod = @views zeros(loop_times, numof_det)
    phi_tod = @views zeros(loop_times, numof_det)
    
    antisun_axis = @views @SVector [1.0, 0.0, 0.0]
    y_axis = @views @SVector [0.0, 1.0, 0.0]
    z_axis = @views @SVector [0.0, 0.0, 1.0]
    spin_axis = @views @SVector [cosd(SS.alpha), 0.0, sind(SS.alpha)]
    
    
    bore_0 = initial_pointings(SS)
    
    #= Generate the quaternion at the initial position of the boresight and detector orientation. =#
    boresight_0 = @views Quaternion(0.0, bore_0)
    detector_vec_0 = @views @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_orientation_0 = @views Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = @views bore_0 × detector_vec_0
    travel_direction_0 = @views Quaternion(0.0, travel_direction_vec_0)

    @views @inbounds @simd for j = eachindex(SS.FP_theta)
        #= 
        Generating the pointhing quaternion of other detectors in the focal plane 
        by adding an angular offset to the boresight based on the SS.FP_theta and SS.FP_phi information. 
        =#
        
        q_theta_in_FP = quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis)
        q_phi_in_FP = quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP
        
        #pointing = @views rotate_by_quarternion!(boresight_0, quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis), quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0))
        
        @views @inbounds @threads for i = eachindex(time_array)
            t = time_array[i]
            #= Generate the quaternion of revolution, precession, and spin at time t. =#
            q_revol = quaternion_rotator(omega_revol, t, z_axis)
            q_prec = quaternion_rotator(omega_prec, t, antisun_axis)
            q_spin = quaternion_rotator(omega_spin, t, spin_axis)
            #= 
            P has the value of the pointing vector (x,y,z) at time t as its imaginary part.
            The direction that the detector orientation is facing is then calculated as det_vec_t.
            =#
            Q = q_revol * q_prec * q_spin
            P = Q * pointing / Q
            q_det = Q * detector_orientation_0 / Q
            travel_direction = Q * travel_direction_0 / Q
            pointing_t = @SVector [P.q1, P.q2, P.q3]
            travel_direction_vec = @SVector [travel_direction.q1, travel_direction.q2, travel_direction.q3]
            #=
            The direction of movement can be calculated by the outer product of the pointhing vector and the detector orientaion vector.
            The vector of meridians can be calculated by combining the outer product of pointing and z axis.
            =#
            longitude = pointing_t × (pointing_t × z_axis)
            ang_t = vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            theta_tod[i, j] = @views ang_t[1]
            phi_tod[i, j] = @views ang_t[2]
            
            #pix_tod[i, j] = ang2pixRing(resol, ang_t[1], ang_t[2])
            #pix_tod[i, j] = ang2pix(_m_, bore_ang[1], bore_ang[2])
            
            #=
            The vector standing perpendicular to the sphere is defined as the outer product of the meridian and the moving_direction. 
            The direction of divergent_vec depends on whether the trajectory enters the sphere from the left or the right side of the meridian.
            =#
            divergent_vec = longitude × travel_direction_vec
            cosK = dot(travel_direction_vec, longitude) / (norm(travel_direction_vec) * norm(longitude))
            #=
            An if statement to prevent errors due to rounding errors that may result in |cosK|>1
            =#
            cosK = ifelse(abs(cosK) > 1.0, sign(cosK), cosK)
            
            psi_tod[i, j] = acos(cosK) * sign(divergent_vec[3]) * sign(pointing_t[3])
        end
    end
    pointings = Dict{String,AbstractArray{Float64}}(
        "theta" => theta_tod,
        "phi" => phi_tod,
        "psi" => psi_tod,
        "time" => time_array
    )
    return pointings#(theta_tod, phi_tod, psi_tod, time_array)
end


@inline function get_pointing_pixels(SS::ScanningStrategy, start, stop)
    resol = @views Resolution(SS.nside)

    omega_spin = rpm2angfreq(SS.spin_rpm)
    omega_prec = rpm2angfreq(SS.prec_rpm)
    omega_revol = (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    
    #= Compute the TOD of the hitpixel and the TOD of the detector orientation at a specified sampling rate from time start to stop. =#
    time_array = @views Vector(start:1/SS.sampling_rate:stop-1/SS.sampling_rate)
    loop_times = length(time_array)
    numof_det = length(SS.FP_theta)
    
    pix_tod = @views zeros(Int64, loop_times, numof_det)
    psi_tod = @views zeros(loop_times, numof_det)
    
    antisun_axis = @views @SVector [1.0, 0.0, 0.0]
    y_axis = @views @SVector [0.0, 1.0, 0.0]
    z_axis = @views @SVector [0.0, 0.0, 1.0]
    spin_axis = @views @SVector [cosd(SS.alpha), 0.0, sind(SS.alpha)]
    
    bore_0 = @views initial_pointings(SS)
    
    #= Generate the quaternion at the initial position of the boresight and detector orientation. =#
    boresight_0 = @views Quaternion(0.0, bore_0)
    detector_vec_0 = @views @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_orientation_0 = @views Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = @views bore_0 × detector_vec_0
    travel_direction_0 = @views Quaternion(0.0, travel_direction_vec_0)

    @inbounds @simd for j = eachindex(SS.FP_theta)
        #= 
        Generating the pointhing quaternion of other detectors in the focal plane 
        by adding an angular offset to the boresight based on the FP_theta and FP_phi information. 
        =#
        q_theta_in_FP = @views quaternion_rotator(deg2rad(SS.FP_theta[j]), 1, y_axis)
        q_phi_in_FP = @views quaternion_rotator(deg2rad(SS.FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP
        @inbounds @threads for i = eachindex(time_array)
            t = @views time_array[i]
            #= Generate the quaternion of revolution, precession, and spin at time t. =#
            q_revol = @views quaternion_rotator(omega_revol, t, z_axis)
            q_prec = @views quaternion_rotator(omega_prec, t, antisun_axis)
            q_spin = @views quaternion_rotator(omega_spin, t, spin_axis)
            #= 
            P has the value of the pointing vector (x,y,z) at time t as its imaginary part.
            The direction that the detector orientation is facing is then calculated as det_vec_t.
            =#
            Q = @views q_revol * q_prec * q_spin
            P = Q * pointing / Q
            q_det = Q * detector_orientation_0 / Q
            travel_direction = Q * travel_direction_0 / Q
            pointing_t = @views @SVector [P.q1, P.q2, P.q3]
            travel_direction_vec = @views @SVector [travel_direction.q1, travel_direction.q2, travel_direction.q3]
            #=
            The direction of movement can be calculated by the outer product of the pointhing vector and the detector orientaion vector.
            The vector of meridians can be calculated by combining the outer product of pointing and z axis.
            =#
            longitude = @views pointing_t × (pointing_t × z_axis)
            ang_t = @views vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            #ang_tod[1, i, j] = ang_t[1]
            #ang_tod[2, i, j] = ang_t[2]
            
            pix_tod[i, j] = @views ang2pixRing(resol, ang_t[1], ang_t[2])
            #pix_tod[i, j] = ang2pix(_m_, bore_ang[1], bore_ang[2])
            
            #=
            The vector standing perpendicular to the sphere is defined as the outer product of the meridian and the moving_direction. 
            The direction of divergent_vec depends on whether the trajectory enters the sphere from the left or the right side of the meridian.
            =#
            divergent_vec = @views longitude × travel_direction_vec
            cosK = dot(travel_direction_vec, longitude) / (norm(travel_direction_vec) * norm(longitude))
            #=
            An if statement to prevent errors due to rounding errors that may result in |cosK|>1
            =#
            cosK = ifelse(abs(cosK) > 1.0, sign(cosK), cosK)
            psi_tod[i, j] = acos(cosK) * sign(divergent_vec[3]) * sign(pointing_t[3])
        end
    end

    return (pix_tod, psi_tod, time_array)
end
