mutable struct scanfield{I<:Integer}
    hitmap::AbstractArray{I,1}
    h::Array
    quantify::DataFrame
    n::Array
    m::Array
    ss::ScanningStrategy
end

function h_nm(field::scanfield, n::Real, m::Real)
    fn(i) = findall(x->x==i, field.n)[1]
    fm(i) = findall(x->x==i, field.m)[1]
    return field.h[fn(n), fm(m),:]
end

function h_nm(hₙₘ::Array, spin_n::Array, spin_m::Array, n::Real, m::Real)
    fn(i) = findall(x->x==i, spin_n)[1]
    fm(i) = findall(x->x==i, spin_m)[1]
    return hₙₘ[fn(n), fm(m),:]
end

function get_hnm_quantify(hₙₘ, spin_n, spin_m)
    df = DataFrame(n            = repeat(spin_n, inner=length(spin_m)),
                   m            = repeat(spin_m, outer=length(spin_n)),
                   mean         = zeros(length(spin_n)*length(spin_m)),
                   std          = zeros(length(spin_n)*length(spin_m)),
                   nanmean      = zeros(length(spin_n)*length(spin_m)),
                   nanstd       = zeros(length(spin_n)*length(spin_m)),
                   nan2one_mean = zeros(length(spin_n)*length(spin_m)),
                   nan2one_std  = zeros(length(spin_n)*length(spin_m)),
    )
    for _n in spin_n
        for _m in spin_m
            row = findfirst((r -> r.n == _n && r.m == _m), eachrow(df))
            abs2_field = abs2.(h_nm(hₙₘ, spin_n, spin_m, _n, _m))
            df[row, :mean]         = mean(abs2_field)
            df[row, :std]          = std(abs2_field)
            df[row, :nanmean]      = nanmean(abs2_field)
            df[row, :nanstd]       = nanstd(abs2_field)
            df[row, :nan2one_mean] = mean(replace(abs2_field, NaN=>1.0))
            df[row, :nan2one_std]  = std(replace(abs2_field, NaN=>1.0))
        end
    end
    return df
end

function get_scanfield(ss::ScanningStrategy,; division::Int, spin_n::Array, spin_m::Array)
    orientation_func_hwp(n, m, ψⱼ, ϕⱼ) = ℯ^(-im*(n*ψⱼ + m*ϕⱼ))
    h      = orientation_func_hwp
    resol  = Resolution(ss.nside)
    npix   = nside2npix(ss.nside)
    chunk  = Int(ss.duration / division)
    ω_hwp  = rpm2angfreq(ss.hwp_rpm)
    hitmap = zeros(Int64, npix)
    hₙₘ     = zeros(Complex{Float64}, (length(spin_n), length(spin_m), npix))
    BEGIN  = 0
    progress = Progress(division)
    @views @inbounds for i = 1:division
        END = i * chunk
        theta, phi, psi, time = get_pointings(ss, BEGIN, END)
        @views @inbounds for j = eachindex(ss.quat)
            theta_j = theta[:,j]
            phi_j   = phi[:,j]
            psi_j   = psi[:,j]
            polang  = get_pol_angle(ss, j)
            @views @inbounds for k = eachindex(time)
                t = time[k]
                p = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t)+polang)
                hitmap[p.Ω] += 1
                @views @inbounds for _n in eachindex(spin_n)
                    @views @inbounds for _m in eachindex(spin_m)
                        hₙₘ[_n, _m, p.Ω] += h(spin_n[_n], spin_m[_m], p.ψ, p.ϕ)
                    end
                end
            end
        end
        BEGIN = END
        next!(progress)
    end
    @views for _n in eachindex(spin_n)
        @views for _m in eachindex(spin_m)
            hₙₘ[_n,_m,:] ./= hitmap
        end
    end
    df = get_hnm_quantify(hₙₘ, spin_n, spin_m)
    return scanfield(hitmap, hₙₘ, df, spin_n, spin_m, ss)
end



function get_scanfield_step(ss::ScanningStrategy,; division::Int, spin_n::Array, spin_m::Array, step)
    orientation_func_hwp(n, m, ψⱼ, ϕⱼ) = ℯ^(-im*(n*ψⱼ + m*ϕⱼ))
    h      = orientation_func_hwp
    resol  = Resolution(ss.nside)
    npix   = nside2npix(ss.nside)
    chunk  = Int(ss.duration / division)
    ω_hwp  = rpm2angfreq(ss.hwp_rpm)
    hitmap = zeros(Int64, npix)
    hₙₘ     = zeros(Complex{Float64}, (length(spin_n), length(spin_m), npix))
    BEGIN  = 0
    progress = Progress(division)
    @views @inbounds for i = 1:division
        END = i * chunk
        theta, phi, psi, time = get_pointings(ss, BEGIN, END, step)
        @views @inbounds for j = eachindex(ss.quat)
            theta_j = theta[:,j]
            phi_j   = phi[:,j]
            psi_j   = psi[:,j]
            polang  = get_pol_angle(ss, j)
            @views @inbounds for k = eachindex(time)
                t = time[k]
                p = pointings(resol, theta_j[k], phi_j[k], psi_j[k], mod2pi(ω_hwp*t)+polang)
                hitmap[p.Ω] += 1
                @views @inbounds for _n in eachindex(spin_n)
                    @views @inbounds for _m in eachindex(spin_m)
                        hₙₘ[_n, _m, p.Ω] += h(spin_n[_n], spin_m[_m], p.ψ, p.ϕ)
                    end
                end
            end
        end
        BEGIN = END
        next!(progress)
    end
    @views for _n in eachindex(spin_n)
        @views for _m in eachindex(spin_m)
            hₙₘ[_n,_m,:] ./= hitmap
        end
    end
    df = get_hnm_quantify(hₙₘ, spin_n, spin_m)
    return scanfield(hitmap, hₙₘ, df, spin_n, spin_m, ss)
end
