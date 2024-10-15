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
    resol = Resolution(ss.nside)
    npix  = nside2npix(ss.nside)
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
        BEGIN = END
        next!(p)
    end
    return (map_div, idx, psi_db, time_db)
end


function get_psi_database(ss::ScanningStrategy,; division::Int, idx, map_div)
    resol = Resolution(ss.nside)
    npix  = nside2npix(ss.nside)
    map_division = map_div*ss.nside
    println("To get a psi database in full-sky, specify idx in the range of [1, $map_division] and run the job.")

    split   = npix/map_division
    month   = Int(ss.duration / division)
    ω_hwp   = rpm2angfreq(ss.hwp_rpm)
    psi_db  = [Float32[] for i in 1:split]
    time_db = [Float32[] for i in 1:split]
    under   = Int(split*(idx-1))
    upper   = Int(split*idx)

    println("ipix_range = [$(under+1),$upper]")
    BEGIN    = 0
    progress = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        theta, phi, psi, time = get_pointings(ss, BEGIN, END)
        @views @inbounds for j = eachindex(ss.quat)
            theta_j   = theta[:,j]
            phi_j     = phi[:,j]
            psi_j     = psi[:,j]
            polang    = get_pol_angle(ss, j)
            @views @inbounds @simd for k = eachindex(time)
                t = time[k]
                p = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t)+polang)
                if under < p.Ω <= upper
                    push!(psi_db[p.Ω - under],  Float32(p.ψ))
                    push!(time_db[p.Ω - under], Float32(t))
                end
            end
        end
        BEGIN = END
        next!(progress)
    end
    return (map_div, idx, psi_db, time_db)
end


w(ψ,ϕ) = @SMatrix [1 cos(2ψ+4ϕ) sin(2ψ+4ϕ)]

function binned_mapmake(ss::ScanningStrategy, division::Int, inputmap::PolarizedHealpixMap, signal)
    resol = Resolution(ss.nside)
    npix = resol.numOfPixels
    chunk = Int(ss.duration / division)
    ω_hwp = rpm2angfreq(ss.hwp_rpm)
    total_signal = zeros(3, 1, npix)
    hitmap = zeros(Int64, npix)
    hitmatrix = zeros(3, 3, npix)
    outmap = zeros(3, 1, npix)
    progress = Progress(division)
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    BEGIN = 0
    @inbounds @views for i = 1:division
        END = i * chunk
        theta, phi, psi, time = get_pointings(ss, BEGIN, END)
        @views @inbounds for j = eachindex(length(ss.quat))
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
    normarize!(total_signal, hitmap)
    normarize!(hitmatrix, hitmap)
    @inbounds @threads for j = eachindex(hitmatrix[1, 1, :])
        det_value = det(hitmatrix[:, :, j])
        if !isnan(det_value) && !isinf(det_value) && det_value != 0
            outmap[:, :, j] = hitmatrix[:, :, j] \ total_signal[:, :, j]
        end
    end
    outmap = transpose([outmap[1,1,:] outmap[2,1,:] outmap[3,1,:]])
    return outmap, hitmap
end
