mutable struct SysPointing{T<:AbstractFloat, AA<:AbstractArray{T}}
    offset_ρ::AA
    offset_χ::AA
    disturbance::AA
end

mutable struct SysGain{T<:AbstractFloat, AA<:AbstractArray{T}}
    offset::AA
    disturbance::AA
end


mutable struct SysEffect
    Gain::SysGain
    Pointing::SysPointing
end

mutable struct InputInfo
    Inputmap::PolarizedHealpixMap
    Systematics::SysEffect
end

function set_input(;
        inputmap,
        gain_offset = [0.0],
        gain_disturbance = [0.0],
        pointing_offset_rho = [0.0],
        pointing_offset_chi = [0.0],
        pointing_disturbance =　[0.0],
    )
    pointing = SysPointing(pointing_offset_rho, pointing_offset_chi, pointing_disturbance)
    gain = SysGain(gain_offset, gain_disturbance)
    systematics = SysEffect(gain, pointing)
    
    inputmap = inputmap |> convert_maps
    InputInfo(inputmap, systematics)
end


@inline function H_μ(t, ω_HWP)
    H = @SMatrix [
        1  0               0
        0  cos(4*ω_HWP*t)  sin(4*ω_HWP*t)
        0  sin(4*ω_HWP*t)  -cos(4*ω_HWP*t)
    ]
    return H
end

@inline function hit_matrix(ρ)
    W = (1.0/4.0).* @SMatrix [
        1.0    cos(ρ)       sin(ρ)
        cos(ρ) cos(ρ)^2     cos(ρ)sin(ρ)
        sin(ρ) cos(ρ)sin(ρ) sin(ρ)^2
    ]
end

function c2d(cl)
    lmax = length(cl)-1
    l = 0:lmax
    lp1 = 1:lmax+1
    return @. (l * lp1)*cl/(2π)
end

w_μ(ψ) = (1.0 / 2.0) .* @SMatrix [1.0; cos(2ψ); sin(2ψ)]

orientation_func(n, ψⱼ) = ℯ^(im*n*ψⱼ)

function hmat(ψⱼ)
    h = orientation_func
    M = @SMatrix [
        1              (1/2)h(2, ψⱼ) (1/2)h(-2, ψⱼ);
        (1/2)h(2,  ψⱼ) (1/4)h(4, ψⱼ) (1/4)         ;
        (1/2)h(-2, ψⱼ) (1/4)         (1/4)h(-4, ψⱼ)
    ]
end


function get_hn_map(SS::ScanningStrategy,; division::Int, nmax=5)
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    
    month = Int(SS.duration / division)
    ω_hwp = rpm2angfreq(SS.hwp_rpm)
    
    hit_map = zeros(npix)
    hₙ = zeros(Complex{Float64}, npix, nmax)
    
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(SS, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp==0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                psi = mod2pi(2*ω_hwp*t) - psi_tod_jth_det[k]
                #psi = psi_tod_jth_det[k]
                hit_map[ipix] += 1
                @views @inbounds for n in 1:nmax
                    hₙ[ipix, n] += cos(n*psi) + sin(n*psi)im
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    @views for n in 1:nmax
        hₙ[:,n] ./= hit_map
    end
    
    outmap = Dict{String, Array}(
        "hitmap" => hit_map,
        "h" => hₙ,
        )
    return outmap
end

#=
function ScanningStrategy2MapInfo(ss::ScanningStrategy, division::Int)
    """
    Nonte:
        This function assumed 2 detectors which is seeing same polarization angle.
    """
    resol = Resolution(ss.nside)
    npix = nside2npix(ss.nside)
    
    month = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    hit_matrix = zeros(npix, 3, 3)
    detmap = zeros(npix)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            #psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                psi = psi_tod_jth_det[k]
                hwp_ang = 4ω_hwp*t
                
                hit_map[ipix] += 2 # Assumed 2 detectors
                Cross[ipix,1,1] += 2sin(hwp_ang - psi) # Assumed 2 detectors
                Cross[ipix,1,2] += 2cos(hwp_ang - psi)
                Cross[ipix,2,1] += 2sin(hwp_ang - 2psi)
                Cross[ipix,2,2] += 2cos(hwp_ang - 2psi)
                Cross[ipix,3,1] += 2sin(hwp_ang - 3psi)
                Cross[ipix,3,2] += 2cos(hwp_ang - 3psi)
                Cross[ipix,4,1] += 2sin(hwp_ang - 4psi)
                Cross[ipix,4,2] += 2cos(hwp_ang - 4psi)
                hit_matrix[ipix,:,:] .+= 2HitMatrix(hwp_ang - 2psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    @views @inbounds @threads for i in eachindex(hit_map)
        Inv_hitmat = inv.(hit_matrix[i,:,:])
        hit_matrix[i,:,:] .= Inv_hitmat
        detmap[i] = det(Inv_hitmat)./hit_map[i]
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    outmap = Dict{String,AbstractArray{Float64}}(
        "hitmap" => hit_map,
        "xl1" => link1,
        "xl2" => link2,
        "xl3" => link3,
        "xl4" => link4,
        "detM" => detmap,
        "dI" => hit_matrix[:,1,1],
        "dQ" => hit_matrix[:,2,2],
        "dU" => hit_matrix[:,3,3]
        )
    return outmap
end
=#

function get_psiDataBase(ss::ScanningStrategy,; division::Int, idx, map_div)
    resol = Resolution(ss.nside)
    npix = nside2npix(ss.nside)
    map_division = map_div*ss.nside
    println("To get a psi database in full-sky, specify idx in the range of [1, $map_division] and run the job.")
    
    split = npix/map_division
    month = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    psi_db = [Float64[] for i in 1:split]
    under = Int(split*(idx-1))
    upper = Int(split*idx)
    println("ipix_range = [$(under+1),$upper]")
    
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k] 
                psi = 2ω_hwp*t - psi_tod_jth_det[k]
                
                if under < ipix <= upper
                    push!(psi_db[ipix-under], psi)
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    return (map_div, idx, psi_db)
end


function get_psi_time_DataBase(ss::ScanningStrategy,; division::Int, idx, map_div)
    #=
    divisionはtodの分割計算
    idxはマップ(npix)をRingOrderでmap_div*nside個の領域に分割したときに北極から順に貼られるインデックス
    map_divはマップをmap_div*nside個に分割するためのパラメータ(4くらいがいい)
    =#
    
    resol = Resolution(ss.nside)
    npix = nside2npix(ss.nside)
    map_division = map_div*ss.nside
    println("To get a psi database in full-sky, specify idx in the range of [1, $map_division] and run the job.")
    
    split = npix/map_division
    month = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    
    psi_db = [Float32[] for i in 1:split]
    time_db = [Float32[] for i in 1:split]
    under = Int(split*(idx-1))
    upper = Int(split*idx)
    println("ipix_range = [$(under+1),$upper]")
    
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(ss, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k] 
                psi = 2ω_hwp*t - psi_tod_jth_det[k]
                
                if under < ipix <= upper
                    push!(psi_db[ipix-under], Float32(psi))
                    push!(time_db[ipix-under], Float32(t))
                end
            end
        end
        #println("Data size psi: ", sizeof(psi_db))
        #println("Data size time: ", sizeof(time_db))
        BEGIN = END
        next!(p)
    end
    return (map_div, idx, psi_db, time_db)
end

w(ψ,ϕ) = @SMatrix [1 cos(2ψ+4ϕ) sin(2ψ+4ϕ)]

function binned_mapmake(ss::ScanningStrategy, division::Int, inputinfo::Falcons.InputInfo, signal)
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
    @inbounds @views for i = 1:division
        END = i * chunk
        theta, phi, psi, time = get_pointings_tuple(ss, BEGIN, END)
        @views @inbounds for j = eachindex(length(ss.FP_theta))
            theta_j = theta[:,j]
            phi_j = phi[:,j]
            psi_j = psi[:,j]
            @inbounds @views for k = eachindex(time)
                t = time[k]
                p = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t))
                dₖ = signal(p, inputmap)
                total_signal[:, :, p.Ω] .+= @SMatrix [dₖ; dₖ*cos(2p.ψ+4p.ϕ); dₖ*sin(2p.ψ+4p.ϕ)]
                hitmap[p.Ω] += 1
                hitmatrix[:, :, p.Ω] .+= transpose(w(p.ψ, p.ϕ)) * w(p.ψ, p.ϕ)
            end
        end
        BEGIN = END
        next!(progress)
    end
    normarize!(resol, total_signal, hitmap)
    normarize!(resol, hitmatrix, hitmap)
    @inbounds @threads for j = eachindex(hitmatrix[1, 1, :])
        det_value = det(hitmatrix[:, :, j])
        if !isnan(det_value) && !isinf(det_value) && det_value != 0
            outmap[:, :, j] = hitmatrix[:, :, j] \ total_signal[:, :, j]
        end
    end
    outmap = transpose([outmap[1,1,:] outmap[2,1,:] outmap[3,1,:]])
    return outmap, hitmap
end

#=
function signal(p::pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    maps.i[p.Ω] + maps.q[p.Ω]*cos(2p.ψ+4p.ϕ) + maps.u[p.Ω]*sin(2p.ψ+4p.ϕ)
end

function interp_signal(p::pointings, maps::PolarizedHealpixMap)
    i = interpolate(maps.i, p.θ, p.φ)
    q = interpolate(maps.q, p.θ, p.φ)
    u = interpolate(maps.u, p.θ, p.φ)
    return i + q*cos(2p.ψ+4p.ϕ) + u*sin(2p.ψ+4p.ϕ)
end
pixbuf = Array{Int}(undef, 4)
weightbuf = Array{Float64}(undef, 4)

function interp_signal(p::pointings, maps::PolarizedHealpixMap, pixbuf, weightbuf)
    i = interpolate(maps.i, p.θ, p.φ, pixbuf, weightbuf)
    q = interpolate(maps.q, p.θ, p.φ, pixbuf, weightbuf)
    u = interpolate(maps.u, p.θ, p.φ, pixbuf, weightbuf)
    return i + q*cos(2p.ψ+4p.ϕ) + u*sin(2p.ψ+4p.ϕ)
end
=#