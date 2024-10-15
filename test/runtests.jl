using Falcons
using Healpix
using JSON
using Test

ss = gen_ScanningStrategy()
ss.sampling_rate = 1.0
ss.nside = 16
ss.duration = 100

function test_get_pointings(;save=false)
    pointings = get_pointings(ss, 0, ss.duration)
    theta = pointings[1][:,1]
    phi = pointings[2][:,1]
    psi = pointings[3][:,1]
    time = pointings[4][:,1]
    if save
        open("reference/pointings.json", "w") do io
            JSON.print(io, Dict("theta" => theta, "phi" => phi, "psi" => psi, "time" => time))
        end
    end
    json = JSON.parsefile("reference/pointings.json")
    @testset "get_pointings-Test" begin
        @test theta == json["theta"]
        @test phi == json["phi"]
        @test psi == json["psi"]
        @test time == json["time"]
    end
end

function test_get_scanfield(;save=false)
    spin_n = [-1, 0, 1/2, 2]
    spin_m = [-1/2, 1/2]
    field = get_scanfield(ss, division=1, spin_n=spin_n, spin_m=spin_m)
    if save
        for n in spin_n
            for m in spin_m
                healpix_map = HealpixMap{Float64, RingOrder}(ss.nside)
                field_non_nan = replace(abs2.(h_nm(field,n,m)), NaN=>0.0)
                healpix_map.pixels .= field_non_nan
                saveToFITS(healpix_map, "!xlink_n=$(n)_m=$(m).fits")
            end
        end
    end
    @testset "get_scanfield-Test" begin
        for n in spin_n
            for m in spin_m
                healpix_map_ref = readMapFromFITS("xlink_n=$(n)_m=$(m).fits", 1, Float64)
                field_non_nan = replace(abs2.(h_nm(field,n,m)), NaN=>0.0)
                @test healpix_map_ref == field_non_nan
            end
        end
    end
end

test_get_pointings(save=false)
test_get_scanfield(save=false)
