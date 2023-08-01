function arcmin2rad(x)
    return deg2rad(x/60.)
end

function pointings(resol::Resolution, θ, φ, ψ, ϕ)
    vec = ang2vec(θ, φ)
    Ω = ang2pixRing(resol, θ, φ)
    ξ = mod2pi(2ϕ) + ψ
    return pointings(vec[1],vec[2],vec[3], θ, φ, Ω, ψ, ϕ, ξ)
end

function true_signal(p::pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    maps.i[p.Ω] + maps.q[p.Ω]*cos(2p.ξ) + maps.u[p.Ω]*sin(2p.ξ)
end

function true_signal(p::pointings, maps::PolarizedHealpixMap)
    maps.i[p.Ω] + maps.q[p.Ω]*cos(2p.ξ) + maps.u[p.Ω]*sin(2p.ξ)
end

function interp_signal(p::pointings, maps::PolarizedHealpixMap)
    i = interpolate(maps.i, p.θ, p.φ)
    q = interpolate(maps.q, p.θ, p.φ)
    u = interpolate(maps.u, p.θ, p.φ)
    return i + q*cos(2p.ξ) + u*sin(2p.ξ)
end

function interp_signal(p::pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    i = interpolate(maps.i, p.θ, p.φ, pixbuf, weightbuf)
    q = interpolate(maps.q, p.θ, p.φ, pixbuf, weightbuf)
    u = interpolate(maps.u, p.θ, p.φ, pixbuf, weightbuf)
    return i + q*cos(2p.ξ) + u*sin(2p.ξ)
end
#=
function normarize!(resol::Resolution, maps::Array, hitmap::Array)
    N = size(maps)[1]
    M = size(maps)[2]
    if N != M
        for i in 1:N
            maps[i,1,:] ./= hitmap
        end
    end
    if N == M
        for i in 1:N
            for j in 1:N
                maps[i,j,:] ./= hitmap
            end
        end
    end
    return maps
end
=#
function gen_signalfield(resol::Resolution, maps::PolarizedHealpixMap)
    alm_i = hp.map2alm(maps.i.pixels)
    alm_q = hp.map2alm(maps.q.pixels)
    alm_u = hp.map2alm(maps.u.pixels)
    di = hp.alm2map_der1(alm_i, resol.nside)
    dq = hp.alm2map_der1(alm_q, resol.nside)
    du = hp.alm2map_der1(alm_u, resol.nside)
    return di,dq,du
end
    
function tayler_expanded_signal(p::pointings, m, II::Falcons.InputInfo)
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

#= multi dets =#
function get_det(path)
    bolonames = []
    open(path) do file
        for (index, line) in enumerate(eachline(file))
            # 行をタブ文字で分割して配列に格納
            row = split(line, '\t')
            filter!(x -> x ≠ "", row)
            # 最初の行をヘッダーとして扱う
            if index == 1
                #ヘッダーをDataFrameの列名として設定
                header = split(line)[2:end]
            else
                push!(bolonames, row[6])
            end
        end
    end
    return bolonames
end


function quat(ϕ, rotate_axis)
    Quaternion([cos(ϕ/2.), rotate_axis[1]*sin(ϕ/2.), rotate_axis[2]*sin(ϕ/2.), rotate_axis[3]*sin(ϕ/2.)])
end

mutable struct OffsetAngle
    #= Offset angle should be defined at apperture coordinate =#
    x::Float64
    y::Float64
    z::Float64
end
#=
function imo2scan_coordinate_offset(ss::ScanningStrategy_imo, offset::OffsetAngle)
    ex = [1., 0., 0.]
    ey = [0., 1., 0.]
    ez = [0., 0., 1.]
    q_boresight = Quaternion(0.,0.,0.,1.)
    q_pol    = Quaternion(0.,1.,0.,0.)
    q_scan   = Quaternion(0.,0.,1.,0.)
    q_d      = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_pol_d  = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    
    spin_axis = [cosd(ss.alpha), 0., sind(ss.alpha)]
    anti_sun_axis = ex
    position = 0.
    
    if ss.start_point == "pole"
        position = π
    end
    for i in eachindex(ss.quat)
        q_imo      = Quaternion(ss.quat[i][4], ss.quat[i][1], ss.quat[i][2], ss.quat[i][3])
        q_offset   = quat(offset.x, ex) * quat(offset.y, ey) * quat(offset.z, ez)
        Q_sun      = quat(position, spin_axis) * quat(deg2rad(90.0-ss.alpha+ss.beta), ey) * q_offset * quat(-π/2., ez) * q_imo
        q_d[i]     = Q_sun * q_boresight / Q_sun
        q_pol_d[i] = Q_sun * q_pol / Q_sun
    end
    
    return (q_d, q_pol_d, q_scan, spin_axis, anti_sun_axis)
end
=#

function imo2scan_coordinate(ss::ScanningStrategy_imo, offset::OffsetAngle)
    ex = [1.0, 0.0, 0.0]
    ey = [0.0, 1.0, 0.0]
    ez = [0.0, 0.0, 1.0]
    ex_apperture_coord = ey
    ey_apperture_coord = ex
    ez_apperture_coord = -ez
    q_boresight = Quaternion(0.,0.,0.,1.)
    q_pol    = Quaternion(0.,1.,0.,0.)
    q_scan   = Quaternion(0.,0.,1.,0.)
    q_d      = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_pol_d  = [Quaternion{Float64}(I) for i in eachindex(ss.quat)]
    q_offset   = quat(offset.x, ex_apperture_coord) * quat(offset.y, ey_apperture_coord) * quat(offset.z, ez_apperture_coord)
    spin_axis = [cosd(ss.alpha), 0, sind(ss.alpha)]
    anti_sun_axis = ex
    q_gamma = 0
    if length(ss.quat) == 1
        if ss.name[1] == "boresight"
            if ss.start_point == "pole"
                flip = π
            elseif ss.start_point == "equator"
                flip = 0
            end
            Q          = quat(deg2rad(90.0-ss.alpha), ey) * quat(flip, ez) * quat(deg2rad(ss.beta), ey) * q_offset
            q_d[1]     = Q * q_boresight / Q
            q_pol_d[1] = Q * q_pol / Q
            return (q_d, q_pol_d, q_scan, spin_axis, anti_sun_axis)
        end
    else
        for i in eachindex(ss.quat)
            telescope  = split.(ss.name[i], "_")[1]
            q_imo      = Quaternion(ss.quat[i][4], ss.quat[i][1], ss.quat[i][2], ss.quat[i][3])
            if telescope     == "000" #LFT
                q_gamma = quat(deg2rad(270), ez)
            elseif telescope == "001" #MFT
                q_gamma = quat(deg2rad(240), ez)
            elseif telescope == "002" #HFT
                q_gamma = quat(deg2rad(30), ez)
            end  
            Q           = quat(deg2rad(ss.beta), ey) *  q_offset *  q_gamma * q_imo
            if telescope == "001" #MFT
                Q = quat(π, ez) * Q
            elseif telescope == "002" #HFT
                Q = quat(π, ez) * Q
            end 
            Q          = quat(deg2rad(90.0-ss.alpha), ey) * Q 
            q_d[i]     = Q * q_boresight / Q
            q_pol_d[i] = Q * q_pol / Q
        end
    end
    return (q_d, q_pol_d, q_scan, spin_axis, anti_sun_axis)
end

function get_pointings_offset(ss::ScanningStrategy_imo, offset::OffsetAngle, start, stop)
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
    
    ey = @SVector [0.0, 1.0, 0.0]
    ez = @SVector [0.0, 0.0, 1.0]
    if ss.coord == "G"
        ey = ecliptic2galactic(ey)
        ez = ecliptic2galactic(ez)
    end

    qb₀, qd₀, qu₀, spin_axis, antisun_axis = imo2scan_coordinate(ss, offset)
    
    @views @inbounds for j = eachindex(ss.quat)
        qp₀ⱼ = qb₀[j]
        @views for i = eachindex(time_array)
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
            
            #θ, ϕ = vec2ang_ver2(p[1], p[2], p[3])
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


function sim_pointing_systematics(ss::ScanningStrategy_imo, 
        division::Int, 
        inputinfo::Falcons.InputInfo,
        offset::OffsetAngle,;
        signal,
        tod_check=false,
        start_time=0,
        end_time=60*60
    )
    w(ψ,ϕ) = @SMatrix [1 cos(2ψ+4ϕ) sin(2ψ+4ϕ)]
    resol = Resolution(ss.nside)
    npix = resol.numOfPixels
    chunk = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    total_signal = zeros(3, 1, npix)
    hitmap = zeros(Int64, npix)
    hitmatrix = zeros(3, 3, npix)
    outmap = zeros(3, 1, npix)
    progress = Progress(division)
    inputmap = inputinfo.Inputmap
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    BEGIN = 0
    no_offset = OffsetAngle(0., 0., 0.)
    
    @inbounds @views for i = 1:division
        END = i * chunk
        theta, phi, psi, time = get_pointings_offset(ss, no_offset, BEGIN, END)
        theta_e, phi_e, psi_e, time_e = get_pointings_offset(ss, offset, BEGIN, END)
        @views @inbounds for j = eachindex(ss.quat)
            theta_j = theta[:,j]
            phi_j = phi[:,j]
            psi_j = psi[:,j]
            theta_e_j = theta_e[:,j]
            phi_e_j = phi_e[:,j]
            psi_e_j = psi_e[:,j]
            @inbounds @views for k = eachindex(time)
                t = time[k]
                p = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t))
                p_err = pointings(resol, theta_e_j[k], phi_e_j[k], psi_e_j[k], mod2pi(ω_hwp*t))
                dₖ = signal(p_err, inputmap)
                total_signal[:, :, p.Ω] .+= @SMatrix [dₖ; dₖ*cos(2p.ξ); dₖ*sin(2p.ξ)]
                hitmatrix[:, :, p.Ω] .+= transpose(w(p.ψ, p.ϕ)) * w(p.ψ, p.ϕ)
                hitmap[p.Ω] += 1
            end
        end
        BEGIN = END
        next!(progress)
    end
    normarize!(resol, total_signal, hitmap)
    normarize!(resol, hitmatrix, hitmap)
    @inbounds @threads for j = 1:npix
        if det(hitmatrix[:, :, j]) != 0
            outmap[:, :, j] = hitmatrix[:, :, j] \ total_signal[:, :, j]
        end
    end
    outmap = transpose([outmap[1,1,:] outmap[2,1,:] outmap[3,1,:]])
    
    if tod_check == true
        theta, phi, psi, time = get_pointings_offset(ss, no_offset, start_time, end_time)
        theta_e, phi_e, psi_e, time_e = get_pointings_offset(ss, offset, start_time, end_time)
        tod_true = zeros(length(time), length(ss.quat))
        tod_true_interp = zeros(length(time), length(ss.quat))
        tod_err_interp = zeros(length(time), length(ss.quat))
        tod_err = zeros(length(time_e), length(ss.quat))
        pol_ang = zeros(length(time), length(ss.quat))
        for j in eachindex(ss.quat)
            for k in eachindex(time)
                t = time[k]
                p = pointings(resol, theta[k,j], phi[k,j], psi[k,j], ω_hwp*t)
                p_err = pointings(resol, theta_e[k,j], phi_e[k,j], psi_e[k,j], ω_hwp*t)
                
                pol_ang[k,j] = p_err.ξ
                #d_true = signal(p, inputmap)
                #d_true_interp = interp_signal(p, inputmap)
                #d_err = signal(p_err, inputmap)
                
                tod_true[k,j] = signal(p, inputmap)
                tod_err[k,j] = signal(p_err, inputmap)
                tod_true_interp[k,j] = interp_signal(p, inputmap)
                tod_err_interp[k,j] = interp_signal(p_err, inputmap)
            end
        end
        data = Dict(
            "outmap" => outmap,
            "hitmap" => hitmap,
            "theta"  => theta,
            "phi" => phi, 
            "psi" => psi,
            "pol_ang" => pol_ang,
            "time" => time,
            "tod_true" => tod_true,
            "tod_true_interp" => tod_true_interp,
            "tod_err" => tod_err,
            "tod_err_interp" => tod_err_interp,
        )
        return data
    end
    return Dict("outmap" => outmap, "hitmap" => hitmap)
end