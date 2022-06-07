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

function xyztod2hitmap(nside::Int, xyz_tod::Array{T}) where {T}
    hit_map = zeros(Int64, nside2npix(nside))
    resol = Resolution(nside)
    for j in eachindex(xyz_tod[1,1,:])
        x_tod_jth_det = @views xyz_tod[1,:,j]
        y_tod_jth_det = @views xyz_tod[2,:,j]
        z_tod_jth_det = @views xyz_tod[3,:,j]
        @inbounds @simd for k = eachindex(xyz_tod[1,:,1])
            ang = @views vec2ang(x_tod_jth_det[k], y_tod_jth_det[k], z_tod_jth_det[k])
            ipix = @views ang2pixRing(resol, ang[1], ang[2])
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

function array2map(map_array::Array)
    nside = npix2nside(length(map_array))
    m = HealpixMap{Float64, RingOrder}(nside)
    m.pixels .= map_array
    return m
end

function convert_maps(healpy_maps)
    PolarizedHealpixMap{Float64,RingOrder}(
        healpy_maps[1,:], 
        healpy_maps[2,:], 
        healpy_maps[3,:],
    )
end


#=
function signal(θ, ϕ, ψ, ipix, t, ω_HWP, sys::Systematics, ;interp=false)
    if interp == false
        power = w_μ(ψ)' * (H_μ(t, ω_HWP) * @SMatrix [sys.Inputmap.i[ipix]; sys.Inputmap.q[ipix]; sys.Inputmap.u[ipix]])
    else
        interp_I = Healpix.interpolate(sys.Inputmap.i, θ, ϕ)
        interp_Q = Healpix.interpolate(sys.Inputmap.q, θ, ϕ)
        interp_U = Healpix.interpolate(sys.Inputmap.u, θ, ϕ,)
        power = w_μ(ψ)' * (H_μ(t, ω_HWP) * @SMatrix [interp_I; interp_Q; interp_U])
    end
    
    if sys.Gain.offset != 0.0
        power .+= sys.Gain.offset
    end
    
    return power
end
=#


