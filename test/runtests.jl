using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
prec1 = 60 * 60 * 3

ss = gen_ScanningStrategy(duration=year)

@testset "ScanningStrategy_structure-Test" begin
    @test typeof(ss) <: ScanningStrategy
    @show fieldnames(ScanningStrategy)
    println("Test ScanningStrategy ==> ", ss)
end

@testset "get_pointings-Test" begin    
    pointings = get_pointings(ss, 0, 30)
    @test typeof(pointings["theta"]) <: Array
    @test typeof(pointings["phi"]) <: Array
    @test typeof(pointings["psi"]) <: Array
    @show pointings["theta"][1:10]
    @show pointings["phi"][1:10]
    @show pointings["psi"][1:10]
end

@testset "ScanningStrategy2map-Test" begin
    n = 6
    outmap = ScanningStrategy2map(ss, n)
    @test length(outmap[1]) == nside2npix(ss.nside)
    @test typeof(outmap[1]) <: Array
    @show outmap[1][1:5]
    @show outmap[2][1:5]
    @show outmap[3][1:5]
    @show outmap[4][1:5]
    @show outmap[5][1:5]
end

