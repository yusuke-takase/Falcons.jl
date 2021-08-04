function ScanningStrategy2MapInfo(SS::ScanningStrategy, division::Int)
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    
    month = Int(SS.duration / division)
    ω_hwp = rpm2angfreq(SS.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    hit_matrix = zeros(npix, 3, 3)
    detmap = zeros(npix)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
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
                Cross[ipix,1,1] += sin(hwp_ang - psi)
                Cross[ipix,1,2] += cos(hwp_ang - psi)
                Cross[ipix,2,1] += sin(hwp_ang - 2psi)
                Cross[ipix,2,2] += cos(hwp_ang - 2psi)
                Cross[ipix,3,1] += sin(hwp_ang - 3psi)
                Cross[ipix,3,2] += cos(hwp_ang - 3psi)
                Cross[ipix,4,1] += sin(hwp_ang - 4psi)
                Cross[ipix,4,2] += cos(hwp_ang - 4psi)
                hit_matrix[ipix,:,:] .+= HitMatrix(hwp_ang - 2psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    #hitmat_ave = @views hit_matrix./hit_map
    @views @inbounds @threads for i = 1:npix
        detmap[i] = det((hit_matrix[i,:,:])./hit_map[i])
        hit_matrix[i,:,:] .= inv.(hit_matrix[i,:,:])
        #detmap[i] = det(inv(hit_matrix[i,:,:]./hit_map[i]))
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
        "dQ" => hit_matrix[:,2,2],
        "dU" => hit_matrix[:,3,3]
        )
    return outmap
end

function get_psiDataBase(SS::ScanningStrategy,; division::Int, idx, map_div)
    #=
    divisionはtodの分割計算
    idxはマップ(npix)をRingOrderでmap_div*nside個の領域に分割したときに北極から順に貼られるインデックス
    map_divはマップをmap_div*nside個に分割するためのパラメータ(4くらいがいい)
    =#
    
    resol = Resolution(SS.nside)
    npix = nside2npix(SS.nside)
    map_division = map_div*SS.nside
    println("To get a psi database in full-sky, specify idx in the range of [1, $map_division] and run the job.")
    
    split = npix/map_division
    month = Int(SS.duration / division)
    ω_hwp = rpm2angfreq(SS.hwp_rpm)
    #hit_map = zeros(npix)
    
    psi_db = [Float64[] for i in 1:split]
    under = Int(split*(idx-1))
    upper = Int(split*idx)
    println("ipix_range = [$(under+1),$upper]")
    
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
                #hwp_ang = 4ω_hwp*t
                #hit_map[ipix] += 1
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

function TwoTelescopes_ScanningStrategy2map(SS1::ScanningStrategy, SS2::ScanningStrategy, division::Int)
    resol = Resolution(SS1.nside)
    npix = nside2npix(SS1.nside)
    
    month = Int(SS1.duration / division)
    ω_hwp = rpm2angfreq(SS1.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    hit_matrix = zeros(npix, 3, 3)
    detmap = zeros(npix)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod1, psi_tod1, time_array1 = get_pointing_pixels(SS1, BEGIN, END)
        pix_tod2, psi_tod2, time_array2 = get_pointing_pixels(SS2, BEGIN, END)
        pix_tod = @views [pix_tod1 pix_tod2]
        psi_tod = @views [psi_tod1 psi_tod2]
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            #psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array1[k]
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
                hit_matrix[ipix,:,:] .+= HitMatrix(hwp_ang - 2psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    
    @views @inbounds @threads for i = 1:npix
        detmap[i] = det((hit_matrix[i,:,:])./hit_map[i])
        hit_matrix[i,:,:] .= inv.(hit_matrix[i,:,:])
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4, detmap]
    return out_map
end


function TwoTelescopes_ScanningStrategy2MapInfo(SS1::ScanningStrategy, SS2::ScanningStrategy, division::Int)
    resol = Resolution(SS1.nside)
    npix = nside2npix(SS1.nside)
    
    month = Int(SS1.duration / division)
    ω_hwp = rpm2angfreq(SS1.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    hit_matrix = zeros(npix, 3, 3)
    detmap = zeros(npix)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod1, psi_tod1, time_array1 = get_pointing_pixels(SS1, BEGIN, END)
        pix_tod2, psi_tod2, time_array2 = get_pointing_pixels(SS2, BEGIN, END)
        pix_tod = @views [pix_tod1 pix_tod2]
        psi_tod = @views [psi_tod1 psi_tod2]
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            #psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array1[k]
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
                hit_matrix[ipix,:,:] .+= HitMatrix(hwp_ang - 2psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    
    @views @inbounds @threads for i = 1:npix
        detmap[i] = det((hit_matrix[i,:,:])./hit_map[i])
        hit_matrix[i,:,:] .= inv.(hit_matrix[i,:,:])
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
        "dQ" => hit_matrix[:,2,2],
        "dU" => hit_matrix[:,3,3]
        )
    return outmap
end


function ThreeTelescopes_ScanningStrategy2map(SS1::ScanningStrategy, SS2::ScanningStrategy, SS3::ScanningStrategy, division::Int)
    resol = Resolution(SS1.nside)
    npix = nside2npix(SS1.nside)
    
    month = Int(SS1.duration / division)
    ω_hwp = rpm2angfreq(SS1.hwp_rpm)
    
    hit_map = zeros(npix)
    Cross = zeros(npix, 4, 2)
    hit_matrix = zeros(npix, 3, 3)
    detmap = zeros(npix)
    BEGIN = 0
    p = Progress(division)
    @views @inbounds for i = 1:division
        END = i * month
        pix_tod1, psi_tod1, time_array1 = get_pointing_pixels(SS1, BEGIN, END)
        pix_tod2, psi_tod2, time_array2 = get_pointing_pixels(SS2, BEGIN, END)
        pix_tod3, psi_tod3, time_array3 = get_pointing_pixels(SS3, BEGIN, END)
        pix_tod = @views [pix_tod1 pix_tod2 pix_tod3]
        psi_tod = @views [psi_tod1 psi_tod2 psi_tod3]
        @views @inbounds for j = eachindex(psi_tod[1,:])
            pix_tod_jth_det = pix_tod[:,j]
            #psi_tod_jth_det = psi_tod[:,j]
            psi_tod_jth_det = ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(psi_tod[:,1])
                t = time_array1[k]
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
                hit_matrix[ipix,:,:] .+= HitMatrix(hwp_ang - 2psi)
            end
        end
        BEGIN = END
        next!(p)
    end
    @views @inbounds @threads for i = 1:npix
        detmap[i] = det((hit_matrix[i,:,:])./hit_map[i])
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4, detmap]
    return out_map
end

function Mapmaking(SS::ScanningStrategy, division::Int)
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
        theta_tod, phi_tod, psi_tod, time_array = get_pointings_tuple(SS, BEGIN, END)
        @inbounds for j = eachindex(psi_tod[1,:])
            theta_tod_jth_det = @views theta_tod[:,j]
            phi_tod_jth_det = @views phi_tod[:,j]
            psi_tod_jth_det = @views ifelse(ω_hwp == 0.0, -psi_tod[:,j], psi_tod[:,j])
            @views @inbounds @simd for k = eachindex(time_array)
                t = @views time_array[k]
                ipix = ang2pixRing(resol, theta_tod_jth_det[k], phi_tod_jth_det[k])
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
