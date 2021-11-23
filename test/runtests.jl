using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
prec1 = 60 * 60 * 3

ss = gen_ScanningStrategy(duration=day)

@testset "ScanningStrategy_structure-Test" begin
    @test typeof(ss) <: ScanningStrategy
    @show fieldnames(ScanningStrategy)
    println("Test ScanningStrategy ==> ", ss)
end
println()
@testset "get_pointings_tuple-Test" begin    
    pointings = get_pointings_tuple(ss, 0, 30)
    @test typeof(pointings[1]) <: Matrix
    @test typeof(pointings[2]) <: Matrix
    @test typeof(pointings[3]) <: Matrix
    @test typeof(pointings[4]) <: StepRangeLen
    #@show pointings[1][1:10]
    #@show pointings[2][1:10]
    #@show pointings[3][1:10]
    #@show pointings[4][1:10]
end
println()
@testset "get_pointing_pixels-Test" begin    
    pointings = get_pointing_pixels(ss, 0, 30)
    @test typeof(pointings[1]) <: Matrix{Int64}
    @test typeof(pointings[2]) <: Matrix{Float64}
    @test typeof(pointings[3]) <: StepRangeLen
    #@show pointings[1][1:10]
    #@show pointings[2][1:10]
    #@show pointings[3][1:10]
end
println()
@testset "get_pointings-Test" begin    
    pointings = get_pointings(ss, 0, 30)
    @test typeof(pointings["theta"]) <: Array
    @test typeof(pointings["phi"]) <: Array
    @test typeof(pointings["psi"]) <: Array
    #@show pointings["theta"][1:10]
    #@show pointings["phi"][1:10]
    #@show pointings["psi"][1:10]
end
println()
@testset "ScanningStrategy2map-Test" begin
    n = 1
    out = ScanningStrategy2map(ss, n)
    @test typeof(out) <: Vector
    #@show outmap[2][1:5]
    #@show outmap[3][1:5]
    #@show outmap[4][1:5]
    #@show outmap[5][1:5]
end
println()
@testset "get_psiDataBase-Test" begin
    DB = get_psiDataBase(ss, division=4, idx=1, map_div=4);
    @test typeof(DB) <: Tuple
end

