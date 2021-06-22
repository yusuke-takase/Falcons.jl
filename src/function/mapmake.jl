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

function Mapmaking(SS::ScanningStrategy, split_num::Int)
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    
    month = Int(SS.duration / split_num)
    ω_hwp = 2π * (SS.hwp_rpm / 60.0)
    
    hit_map = zeros(npix)
    Cross = zeros(2,4, npix)
    BEGIN = 0
    p = Progress(split_num)
    @views @inbounds for i = 1:split_num
        END = i * month
        pointings = get_pointings(SS, BEGIN, END)
        @inbounds for j = eachindex(pointings["psi"][1,:])
            theta_tod_jth_det = @views pointings["theta"][:,j]
            phi_tod_jth_det = @views pointings["phi"][:,j]
            psi_tod_jth_det = @views ifelse(ω_hwp == 0.0, -pointings["psi"][:,j], pointings["psi"][:,j])
            @views @inbounds @simd for k = eachindex(pointings["psi"][:,1])
                t = @views pointings["time"][k]
                ipix = ang2pixRing(resol, theta_tod_jth_det[k], phi_tod_jth_det[k])
                psi = psi_tod_jth_det[k]
                hwp_ang = 4ω_hwp*t
                
                hit_map[ipix] += 1
                Cross[1,1,ipix] += sin(hwp_ang - psi)
                Cross[2,1,ipix] += cos(hwp_ang - psi)
                Cross[1,2,ipix] += sin(hwp_ang - 2psi)
                Cross[2,2,ipix] += cos(hwp_ang - 2psi)
                Cross[1,3,ipix] += sin(hwp_ang - 3psi)
                Cross[2,3,ipix] += cos(hwp_ang - 3psi)
                Cross[1,4,ipix] += sin(hwp_ang - 4psi)
                Cross[2,4,ipix] += cos(hwp_ang - 4psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    link1 = @views @. (Cross[1,1,:]/hit_map)^2 + (Cross[2,1,:]/hit_map)^2
    link2 = @views @. (Cross[1,2,:]/hit_map)^2 + (Cross[2,2,:]/hit_map)^2
    link3 = @views @. (Cross[1,3,:]/hit_map)^2 + (Cross[2,3,:]/hit_map)^2
    link4 = @views @. (Cross[1,4,:]/hit_map)^2 + (Cross[2,4,:]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4]
    return out_map
end

function ScanningStrategy2map(SS::ScanningStrategy, split_num::Int)
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    
    month = Int(SS.duration / split_num)
    ω_hwp = 2π * (SS.hwp_rpm / 60.0)
    
    hit_map = zeros(npix)
    Cross = zeros((2,4, npix))
    BEGIN = 0
    p = Progress(split_num)
    @views @inbounds for i = 1:split_num
        END = i * month
        pix_tod, psi_tod, time_array = get_pointing_pixels(SS, BEGIN, END)
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            #psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array[k]
                ipix = pix_tod_jth_det[k]
                psi = psi_tod_jth_det[k]
                hwp_ang = 4ω_hwp*t
                
                hit_map[ipix] += 1
                Cross[1,1,ipix] += sin(hwp_ang - psi)
                Cross[2,1,ipix] += cos(hwp_ang - psi)
                Cross[1,2,ipix] += sin(hwp_ang - 2psi)
                Cross[2,2,ipix] += cos(hwp_ang - 2psi)
                Cross[1,3,ipix] += sin(hwp_ang - 3psi)
                Cross[2,3,ipix] += cos(hwp_ang - 3psi)
                Cross[1,4,ipix] += sin(hwp_ang - 4psi)
                Cross[2,4,ipix] += cos(hwp_ang - 4psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    link1 = @views @. (Cross[1,1,:]/hit_map)^2 + (Cross[2,1,:]/hit_map)^2
    link2 = @views @. (Cross[1,2,:]/hit_map)^2 + (Cross[2,2,:]/hit_map)^2
    link3 = @views @. (Cross[1,3,:]/hit_map)^2 + (Cross[2,3,:]/hit_map)^2
    link4 = @views @. (Cross[1,4,:]/hit_map)^2 + (Cross[2,4,:]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4]
    return out_map
end

function Genmap(map_array::Array)
    nside = npix2nside(length(map_array))
    m = Map{Float64, RingOrder}(nside)
    m.pixels .= map_array
    return m
end