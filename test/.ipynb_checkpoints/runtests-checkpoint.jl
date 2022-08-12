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
    #@test typeof(pointings[4]) <: StepRangeLen
    #@show pointings[1][1:10]
    #@show pointings[2][1:10]
    #@show pointings[3][1:10]
    #@show pointings[4][1:10]
end