using StaticArrays

mutable struct SysPointing{AF<:AbstractFloat}
    offset_ρ::AF
    offset_χ::AF
    disturbance::AF
end

mutable struct SysGain{AF<:AbstractFloat}
    offset::AF
    disturbance::AF
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
        gain_offset = 0.0,
        gain_disturbance = 0.0,
        pointing_offset_rho = 0.0,
        pointing_offset_chi = 0.0,
        pointing_disturbance =0.0,
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

function get_hmap(SS::ScanningStrategy,; division::Int, nmax=5)
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


