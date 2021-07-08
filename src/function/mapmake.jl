function angtod2hitmap(nside::Int, theta_tod::Array{T}, phi_tod::Array{T}) where {T}
    hit_map = zeros(Int64, nside2npix(nside))
    resol = Resolution(nside)
    for j in eachindex(theta_tod[1,:])
        theta_tod_jth_det = @views theta_tod[:,j]
        phi_tod_jth_det = @views phi_tod[:,j]
        @inbounds @simd for k = eachindex(theta_tod[:,1])
            ipix = @views ang2pixRing(resol, theta_tod_jth_det[k], phi_tod_jth_det[k])
            hit_map[ipix] += 1
        end
    end
    return hit_map
end

function pixtod2hitmap(nside::Int, pixtod::Array{T}) where {T}
    npix = nside2npix(nside)
    hit_map = zeros(Int64, npix)
    for i = eachindex(pixtod)
        hit_map[pixtod[i]] += 1
    end
    return hit_map
end

function ScanningStrategy2map(SS::ScanningStrategy, division::Int)
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    
    month = Int(SS.duration / division)
    ω_hwp = rpm2angfreq(SS.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(SS, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                psi = psi_tod_jth_det[k]
                hwp_ang = 4ω_hwp*t
                
                hit_map[ipix] += 1
                Cross[ipix,1,1] += sin(hwp_ang - psi)
                Cross[ipix,1,2] += cos(hwp_ang - psi)
                Cross[ipix,2,1] += sin(hwp_ang - 2psi)
                Cross[ipix,2,2] += cos(hwp_ang - 2psi)
                Cross[ipix,3,1] += sin(hwp_ang - 3psi)
                Cross[ipix,3,2] += cos(hwp_ang - 3psi)
                Cross[ipix,4,1] += sin(hwp_ang - 4psi)
                Cross[ipix,4,2] += cos(hwp_ang - 4psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4]
    return out_map
end

function Genmap(map_array::Array)
    nside = npix2nside(length(map_array))
    m = Map{Float64, RingOrder}(nside)
    m.pixels .= map_array
    return m
end

@inline function Hᵢ(t)
    ω_HWP = @views 2 * π * (46 / 60)
    s = @views sin(4*ω_HWP*t)
    c = @views cos(4*ω_HWP*t)
    H = @views @SMatrix [
        1 0  0
        0 c  s
        0 s -c
    ]
    return H
end

@inline function HitMatrix(rho)
    s = sin(rho)
    c = cos(rho)
    D = @views (1/1).* @SMatrix [
        1 c   s
        c c^2 s*c
        s s*c s^2
    ]
    return D
end

@inline function pᵢ(pix_i, psi, t, obsmap)
    _wᵢ(psi) = @SMatrix [1.0/2.0; cos(2psi) / 2.0; sin(2psi) / 2.0]
    I = @views obsmap[1,:]
    Q = @views obsmap[2,:]
    U = @views obsmap[3,:]
    return _wᵢ(psi)' * (Hᵢ(t) * [I[pix_i]; Q[pix_i]; U[pix_i]])
end

_wᵢ(psi) = @SMatrix [1.0/2.0; cos(2psi) / 2.0; sin(2psi) / 2.0]
