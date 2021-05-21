function HitMap(nside, theta_tod, phi_tod)
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


function Mapmaking(scan_strategy_struct, split_num)
    resol = Resolution(scan_strategy_struct.nside)
    npix = nside2npix(scan_strategy_struct.nside)
    
    month = Int(scan_strategy_struct.times / split_num)
    hwp_revol_rate = 2.0 * π * (scan_strategy_struct.hwp_rpm/ 60.0)
    
    hit_map = zeros(Int32, npix)
    Cross = zeros(Float32, (2,4, npix))
    BEGIN = 0
    
    @inbounds @simd for i = 1:split_num
        println("process=", i, "/", split_num)
        END = i * month
        
        theta_tod, phi_tod, psi_tod = get_scan_tod(scan_strategy_struct, BEGIN, END)
        println("Start mapmaking!")
        
        @inbounds @simd for j = eachindex(theta_tod[1,:])
            theta_tod_jth_det = @views theta_tod[:,j]
            phi_tod_jth_det = @views phi_tod[:,j]
            psi_tod_jth_det = @views psi_tod[:,j]
            @inbounds @simd for k = eachindex(theta_tod[:,1])
                ipix = @views ang2pixRing(resol, theta_tod_jth_det[k], phi_tod_jth_det[k])
                psi = @views psi_tod_jth_det[k]
                TIME = @views BEGIN + (k - 1) / scan_strategy_struct.sampling_rate

                hit_map[ipix] += 1

                Cross[1,1,ipix] += sin(psi)
                Cross[2,1,ipix] += cos(psi)
                Cross[1,2,ipix] += sin(2psi)
                Cross[2,2,ipix] += cos(2psi)
                Cross[1,3,ipix] += sin(3psi)
                Cross[2,3,ipix] += cos(3psi)
                Cross[1,4,ipix] += sin(4psi)
                Cross[2,4,ipix] += cos(4psi)

            end
        end
        BEGIN = END + 1
    end
    
    link1 = @views @. (Cross[1,1,:]/hit_map)^2 + (Cross[2,1,:]/hit_map)^2
    link2 = @views @. (Cross[1,2,:]/hit_map)^2 + (Cross[2,2,:]/hit_map)^2
    link3 = @views @. (Cross[1,3,:]/hit_map)^2 + (Cross[2,3,:]/hit_map)^2
    link4 = @views @. (Cross[1,4,:]/hit_map)^2 + (Cross[2,4,:]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4]
    
    return out_map
end
