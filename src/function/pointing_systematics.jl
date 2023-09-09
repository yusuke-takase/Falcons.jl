function true_signal(p::Pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    maps.i[p.Ω] + maps.q[p.Ω]*cos(2p.ξ) + maps.u[p.Ω]*sin(2p.ξ)
end

function true_signal(p::Pointings, maps::PolarizedHealpixMap)
    maps.i[p.Ω] + maps.q[p.Ω]*cos(2p.ξ) + maps.u[p.Ω]*sin(2p.ξ)
end

function interp_signal(p::Pointings, maps::PolarizedHealpixMap)
    i = interpolate(maps.i, p.θ, p.φ)
    q = interpolate(maps.q, p.θ, p.φ)
    u = interpolate(maps.u, p.θ, p.φ)
    return i + q*cos(2p.ξ) + u*sin(2p.ξ)
end


function interp_signal(p::Pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    i = interpolate(maps.i, p.θ, p.φ, pixbuf, weightbuf)
    q = interpolate(maps.q, p.θ, p.φ, pixbuf, weightbuf)
    u = interpolate(maps.u, p.θ, p.φ, pixbuf, weightbuf)
    return i + q*cos(2p.ξ) + u*sin(2p.ξ)
end

function gen_signalfield(resol::Resolution, maps::PolarizedHealpixMap)
    alm_i = hp.map2alm(maps.i.pixels)
    alm_q = hp.map2alm(maps.q.pixels)
    alm_u = hp.map2alm(maps.u.pixels)
    di = hp.alm2map_der1(alm_i, resol.nside)
    dq = hp.alm2map_der1(alm_q, resol.nside)
    du = hp.alm2map_der1(alm_u, resol.nside)
    return di,dq,du
end
    
function tayler_expanded_signal(p::Pointings, m, II::Falcons.InputInfo)
    ∂I = @views m[1][3,p.Ω] - m[1][2,p.Ω]im
    ∂P = @views m[2][3,p.Ω] + m[3][2,p.Ω] - (m[2][2,p.Ω] - m[3][3,p.Ω])im
    ∂̄P = @views m[2][3,p.Ω] - m[3][2,p.Ω] + (m[2][2,p.Ω] + m[3][3,p.Ω])im
    I = II.Inputmap.i.pixels[p.Ω]
    P = II.Inputmap.q.pixels[p.Ω] + II.Inputmap.u.pixels[p.Ω]*im
    ρ = II.Systematics.Pointing.offset_ρ[1]
    χ = II.Systematics.Pointing.offset_χ[1]
    I1 = I - (ρ/2)*(ℯ^(im*(p.ψ+χ))*∂I + ℯ^(-im*(p.ψ+χ))*conj(∂I))
    P1 = (1/2) * (P*ℯ^(im*(2p.ψ+4p.ϕ))        - (ρ/2) * (ℯ^(im*(3p.ψ+4p.ϕ+χ))*∂P       + ℯ^(im*(p.ψ+4p.ϕ-χ))*∂̄P ))
    P2 = (1/2) * (conj(P)*ℯ^(-im*(2p.ψ+4p.ϕ)) - (ρ/2) * (ℯ^(im*(-p.ψ-4p.ϕ+χ))*conj(∂̄P) + ℯ^(im*(-3p.ψ-4p.ϕ-χ))*conj(∂P) ))
    return I1 + P1 + P2
end

mutable struct OffsetAngles{AF<:AbstractFloat}
    #= Offset angle should be defined at apperture coordinate =#
    x::AF
    y::AF
    z::AF
end

function imo2ecl_coordinates(ss::ScanningStrategy, offset::OffsetAngles)
    #=
    This function is loading imo directory. The imo v2.0< have a bug that LFT orientation. 
    In order to modify the bug, the comment of `g_gamma` have to be erase. 
    =#
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    ex_apperture_coord = ey
    ey_apperture_coord = ex
    ez_apperture_coord = -ez
    
    q_boresight = Quaternion(0.,0.,0.,1.)
    q_pol       = Quaternion(0.,1.,0.,0.)
    q_dets      = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_pol_dets  = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_offset    = rotate_quat(offset.x, ex_apperture_coord) * rotate_quat(offset.y, ey_apperture_coord) * rotate_quat(offset.z, ez_apperture_coord)
    q_gamma     = 0.
    polang      = 0.
    for i in eachindex(ss.quat)
        if ss.name[i] == "boresight"
            if ss.start_point == "pole"
                flip = π
            elseif ss.start_point == "equator"
                flip = 0.
            end
            Q              = rotate_quat(deg2rad(90.0-ss.alpha), ey) * rotate_quat(flip, ez) * rotate_quat(deg2rad(ss.beta), ey) * q_offset
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
            Q              = rotate_quat(deg2rad(ss.beta), ey) *  q_offset * q_gamma * q_imo
            q_dets[i]      = Q * q_boresight / Q
            q_pol_dets[i]  = Q * q_pol / Q
        end
    end
    return (q_dets, q_pol_dets)
end

function get_pointings_offset(ss::ScanningStrategy, offset::OffsetAngles, start, stop)
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
    q_scan_direction             = Quaternion(0.,0.,1.,0.)
    q_point, q_pol_angle = imo2ecl_coordinates(ss, offset) # This line is only different from get_pointings()
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
            
            θ, ϕ           = vec2ang_lbsim_phi(vec_point[1], vec_point[2], vec_point[3])
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



function sim_pointing_systematics(ss::ScanningStrategy,
        nside_out::Int,
        division::Int, 
        inputmap::PolarizedHealpixMap,
        offset::OffsetAngles,;
        signal,
        tod_check=false,
        start_time=0,
        end_time=60*60
    )
    w(ψ,ϕ)       = @SMatrix [1 cos(2ψ+4ϕ) sin(2ψ+4ϕ)]
    resol        = Resolution(ss.nside)
    resol_out    = Resolution(nside_out)
    npix         = nside2npix(nside_out) #resol.numOfPixels
    chunk        = Int(ss.duration / division)
    ω_hwp        = rpm2angfreq(ss.hwp_rpm)
    total_signal = zeros(3, 1, npix)
    hitmap       = zeros(Int64, npix)
    hitmatrix    = zeros(3, 3, npix)
    outmap       = zeros(3, 1, npix)
    progress     = Progress(division)
    pixbuf       = Array{Int}(undef, 4)
    weightbuf    = Array{Float64}(undef, 4)
    BEGIN        = 0
    no_offset    = OffsetAngles(0., 0., 0.)
    @inbounds @views for i = 1:division
        END = i * chunk
        theta, phi, psi, time         = get_pointings_offset(ss, no_offset, BEGIN, END)
        theta_e, phi_e, psi_e, time_e = get_pointings_offset(ss, offset, BEGIN, END)
        @views @inbounds for j = eachindex(ss.quat)
            theta_j   = theta[:,j]
            phi_j     = phi[:,j]
            psi_j     = psi[:,j]
            theta_e_j = theta_e[:,j]
            phi_e_j   = phi_e[:,j]
            psi_e_j   = psi_e[:,j]
            polang    = get_pol_angle(ss, j)
            @inbounds @views for k = eachindex(time)
                t     = time[k]
                p     = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t)+polang)
                p_out = pointings(resol_out, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t)+polang)
                p_err = pointings(resol, theta_e_j[k], phi_e_j[k], psi_e_j[k], mod2pi(ω_hwp*t)+polang+offset.z)
                dₖ    = signal(p_err, inputmap)
                total_signal[:, :, p_out.Ω] .+= @SMatrix [dₖ; dₖ*cos(2p.ξ); dₖ*sin(2p.ξ)]
                hitmatrix[:, :, p_out.Ω]    .+= transpose(w(p.ψ, p.ϕ)) * w(p.ψ, p.ϕ)
                hitmap[p_out.Ω] += 1
            end
        end
        BEGIN = END
        next!(progress)
    end
    normarize!(resol, total_signal, hitmap)
    normarize!(resol, hitmatrix, hitmap)
    @inbounds @threads for j = eachindex(hitmatrix[1, 1, :])
        det_value = det(hitmatrix[:, :, j])
        if !isnan(det_value) && det_value != 0
            outmap[:, :, j] = hitmatrix[:, :, j] \ total_signal[:, :, j]
        end
    end
    outmap = transpose([outmap[1,1,:] outmap[2,1,:] outmap[3,1,:]])
    if tod_check == true
        theta, phi, psi, time         = get_pointings_offset(ss, no_offset, start_time, end_time)
        theta_e, phi_e, psi_e, time_e = get_pointings_offset(ss, offset, start_time, end_time)
        tod_true        = zeros(length(time), length(ss.quat))
        tod_true_interp = zeros(length(time), length(ss.quat))
        tod_err_interp  = zeros(length(time), length(ss.quat))
        tod_err         = zeros(length(time), length(ss.quat))
        pol_ang         = zeros(length(time), length(ss.quat))
        pol_ang_err     = zeros(length(time), length(ss.quat))
        hwp_ang         = zeros(length(time), length(ss.quat))
        for j in eachindex(ss.quat)
            for k in eachindex(time)
                t                    = time[k]
                p                    = pointings(resol, theta[k,j], phi[k,j], psi[k,j], ω_hwp*t)
                p_err                = pointings(resol, theta_e[k,j], phi_e[k,j], psi_e[k,j], ω_hwp*t)
                pol_ang[k,j]         = p.ξ
                pol_ang_err[k,j]     = p_err.ξ
                hwp_ang[k,j]         = mod2pi(p.ϕ)
                tod_true[k,j]        = signal(p, inputmap)
                tod_err[k,j]         = signal(p_err, inputmap)
                tod_true_interp[k,j] = interp_signal(p, inputmap)
                tod_err_interp[k,j]  = interp_signal(p_err, inputmap)
            end
        end
        data = Dict(
            "outmap"          => outmap,
            "hitmap"          => hitmap,
            "theta"           => theta,
            "phi"             => phi, 
            "psi"             => psi,
            "pol_ang"         => pol_ang,
            "pol_ang_err"     => pol_ang_err,
            "hwp_ang"         => hwp_ang,
            "time"            => time,
            "tod_true"        => tod_true,
            "tod_true_interp" => tod_true_interp,
            "tod_err"         => tod_err,
            "tod_err_interp"  => tod_err_interp,
        )
        return data
    end
    return Dict("outmap" => outmap, "hitmap" => hitmap)
end