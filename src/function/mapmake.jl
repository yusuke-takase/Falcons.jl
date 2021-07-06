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

function ScanningStrategy2map(SS::ScanningStrategy, division::Int)
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
        #detmap[i] = det(inv(hit_matrix[i,:,:]./hit_map[i]))
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4, detmap]
    return out_map
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
                hwp_ang = 4ω_hwp*t
                #hit_map[ipix] += 1
                if under < ipix <= upper
                    push!(psi_db[ipix-under], psi)
                end
            end
        end
        BEGIN = END
        next!(p)
    end
    return psi_db
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
    end
    link1 = @views @. (Cross[:,1,1]/hit_map)^2 + (Cross[:,1,2]/hit_map)^2
    link2 = @views @. (Cross[:,2,1]/hit_map)^2 + (Cross[:,2,2]/hit_map)^2
    link3 = @views @. (Cross[:,3,1]/hit_map)^2 + (Cross[:,3,2]/hit_map)^2
    link4 = @views @. (Cross[:,4,1]/hit_map)^2 + (Cross[:,4,2]/hit_map)^2
    out_map = @views [hit_map, link1, link2, link3, link4, detmap]
    
    return out_map
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
