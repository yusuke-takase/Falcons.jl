function HitMap(nside, pix)
    npix = nside2npix(nside)
    hit_map = zeros(Int64, npix)
    for i = eachindex(pix)
        hit_map[pix[i]] += 1
    end
    return hit_map
end


function Mapmaking(scan_strategy_struct)
    npix = nside2npix(scan_strategy_struct.nside)
    split_num = 2
    month = Int(scan_strategy_struct.times / split_num)
    hwp_revol_rate = 2.0 * π * (scan_strategy_struct.hwp_rpm/ 60.0)
    
    hit_map = zeros(Int32, npix)
    Cross = zeros(Float32, (2,4, npix))
    BEGIN = 0
    
    @inbounds @simd for i = 1:split_num
        println("process=", i, "/", split_num)
        END = i * month
        
        pix_tod, psi_tod = get_scan_tod(scan_strategy_struct, BEGIN, END)
        println("Start mapmaking!")
        
        @inbounds @simd for k = eachindex(pix_tod)
            bore = @views pix_tod[k]
            psi = @views psi_tod[k]
            TIME = @views BEGIN + (k - 1) / scan_strategy_struct.sampling_rate
            
            hit_map[bore] += 1
            
            Cross[1,1,bore] += sin(psi)
            Cross[2,1,bore] += cos(psi)
            Cross[1,2,bore] += sin(2psi)
            Cross[2,2,bore] += cos(2psi)
            Cross[1,3,bore] += sin(3psi)
            Cross[2,3,bore] += cos(3psi)
            Cross[1,4,bore] += sin(4psi)
            Cross[2,4,bore] += cos(4psi)
            
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
