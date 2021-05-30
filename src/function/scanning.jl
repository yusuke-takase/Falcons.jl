mutable struct ScanningStrategy
    nside::Int
    times::Int
    sampling_rate::Int
    alpha::AbstractFloat
    beta::AbstractFloat
    prec_period::AbstractFloat
    spin_rpm::AbstractFloat
    hwp_rpm::AbstractFloat
    FP_theta:: AbstractArray{AbstractFloat,1}
    FP_phi:: AbstractArray{AbstractFloat,1}
    start_point::AbstractString
    ScanningStrategy() = new()
end


@inline function get_scan_tod(ScanningStrategyStructure, start, stop)
    SSS = @views ScanningStrategyStructure
    nside = SSS.nside
    alpha = @views deg2rad(SSS.alpha)
    beta = @views deg2rad(SSS.beta)
    FP_theta = SSS.FP_theta
    FP_phi = SSS.FP_phi
    omega_spin = @views (2π / 60) * SSS.spin_rpm
    omega_prec = @views (2π / 60) / SSS.prec_period
    
    #= Compute the TOD of the hitpixel and the TOD of the detector orientation at a specified sampling rate from time start to stop. =#
    
    smp_rate = SSS.sampling_rate
    loop_times = @views ((stop - start) * smp_rate)  + smp_rate
    #pix_tod = @views zeros(Int32, loop_times, length(FP_theta))
    psi_tod = @views zeros(Float32, loop_times, length(FP_theta))
    ang_tod = @views zeros(Float32, 2, loop_times, length(FP_theta))
    
    resol = @views Resolution(nside)
    
    antisun_axis = @views @SVector [1.0, 0.0, 0.0]
    spin_axis = @views @SVector [cos(alpha), 0.0, sin(alpha)]
    y_axis = @views @SVector [0.0, 1.0, 0.0]
    z_axis = @views @SVector [0.0, 0.0, 1.0]
    omega_revol = @views (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    
    if SSS.start_point =="pole"
        #= Set the initial position of the boresight to near the zenith. =#
        bore_0 = @views @SVector [sin(alpha-beta), 0, cos(alpha-beta)]
    end
    if SSS.start_point == "equator"
        #= Set the initial position of the boresight to near the equator. =#
        bore_0 = @views @SVector [cos(alpha-beta), 0 , sin(alpha-beta)]
    end
    #= Generate the quaternion at the initial position of the boresight and detector orientation. =#
    boresight_0 = @views Quaternion(0.0, bore_0)
    detector_vec_0 = @views @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_orientation_0 = @views Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = @views bore_0 × detector_vec_0
    travel_direction_0 = @views Quaternion(0.0, travel_direction_vec_0)

    @views @inbounds @simd for j in eachindex(FP_theta)
        #= 
        Generating the pointhing quaternion of other detectors in the focal plane 
        by adding an angular offset to the boresight based on the FP_theta and FP_phi information. 
        =#
        q_theta_in_FP = quaternion_rotator(deg2rad(FP_theta[j]), 1, y_axis)
        q_phi_in_FP = quaternion_rotator(deg2rad(FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP
        @views @inbounds @threads for i = 1:loop_times
            t = start + (i - 1) / smp_rate
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
            bore_ang = vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            ang_tod[1, i, j] = bore_ang[1]
            ang_tod[2, i, j] = bore_ang[2]
            
            #pix_tod[i, j] = ang2pixRing(resol, bore_ang[1], bore_ang[2])
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
            if abs(cosK) > 1.0
                cosK = sign(cosK)
            end
            psi_tod[i, j] = acos(cosK) * sign(divergent_vec[3]) * sign(pointing_t[3])
        end
    end
    return ang_tod[1,:,:], ang_tod[2,:,:], psi_tod
end

@inline function get_scan_tod_pix(ScanningStrategyStructure, start, stop)
    SSS = @views ScanningStrategyStructure
    nside = SSS.nside
    alpha = @views deg2rad(SSS.alpha)
    beta = @views deg2rad(SSS.beta)
    FP_theta = SSS.FP_theta
    FP_phi = SSS.FP_phi
    omega_spin = @views (2π / 60) * SSS.spin_rpm
    omega_prec = @views (2π / 60) / SSS.prec_period
    
    #= Compute the TOD of the hitpixel and the TOD of the detector orientation at a specified sampling rate from time start to stop. =#
    
    smp_rate = SSS.sampling_rate
    loop_times = @views ((stop - start) * smp_rate)  + smp_rate
    pix_tod = @views zeros(Int32, loop_times, length(FP_theta))
    psi_tod = @views zeros(Float32, loop_times, length(FP_theta))
    #ang_tod = @views zeros(Float32, 2, loop_times, length(FP_theta))
    
    resol = @views Resolution(nside)
    
    antisun_axis = @views @SVector [1.0, 0.0, 0.0]
    spin_axis = @views @SVector [cos(alpha), 0.0, sin(alpha)]
    y_axis = @views @SVector [0.0, 1.0, 0.0]
    z_axis = @views @SVector [0.0, 0.0, 1.0]
    omega_revol = @views (2π) / (60.0 * 60.0 * 24.0 * 365.25)
    
    if SSS.start_point =="pole"
        #= Set the initial position of the boresight to near the zenith. =#
        bore_0 = @views @SVector [sin(alpha-beta), 0, cos(alpha-beta)]
    end
    if SSS.start_point == "equator"
        #= Set the initial position of the boresight to near the equator. =#
        bore_0 = @views @SVector [cos(alpha-beta), 0 , sin(alpha-beta)]
    end
    #= Generate the quaternion at the initial position of the boresight and detector orientation. =#
    boresight_0 = @views Quaternion(0.0, bore_0)
    detector_vec_0 = @views @SVector [-boresight_0.q3, 0.0, boresight_0.q1]
    detector_orientation_0 = @views Quaternion(0.0, detector_vec_0)
    travel_direction_vec_0 = @views bore_0 × detector_vec_0
    travel_direction_0 = @views Quaternion(0.0, travel_direction_vec_0)

    @views @inbounds @simd for j in eachindex(FP_theta)
        #= 
        Generating the pointhing quaternion of other detectors in the focal plane 
        by adding an angular offset to the boresight based on the FP_theta and FP_phi information. 
        =#
        q_theta_in_FP = quaternion_rotator(deg2rad(FP_theta[j]), 1, y_axis)
        q_phi_in_FP = quaternion_rotator(deg2rad(FP_phi[j]), 1, bore_0)
        q_for_FP = q_phi_in_FP * q_theta_in_FP
        pointing = q_for_FP * boresight_0 / q_for_FP
        @views @inbounds @threads for i = 1:loop_times
            t = start + (i - 1) / smp_rate
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
            bore_ang = vec2ang_ver2(pointing_t[1], pointing_t[2], pointing_t[3])

            #ang_tod[1, i, j] = bore_ang[1]
            #ang_tod[2, i, j] = bore_ang[2]
            
            pix_tod[i, j] = ang2pixRing(resol, bore_ang[1], bore_ang[2])
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
            if abs(cosK) > 1.0
                cosK = sign(cosK)
            end
            psi_tod[i, j] = acos(cosK) * sign(divergent_vec[3]) * sign(pointing_t[3])
        end
    end
    return pix_tod, psi_tod
end
