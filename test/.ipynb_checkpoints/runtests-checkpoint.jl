using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
prec1 = 60 * 60 * 3

ss = ScanningStrategy()
ss.nside = 128
ss.times = year
ss.sampling_rate = 2
ss.FP_theta = [0.0]
ss.FP_phi = [0.0]
ss.alpha = 55
ss.beta = 60
ss.prec_period = 180.223
ss.spin_rpm = 0.04
ss.hwp_rpm = 0.05
ss.start_point = "pole"

@testset "ScanningStrategy_structure-Test" begin
    @test typeof(ss) <: ScanningStrategy
    @show fieldnames(ScanningStrategy)
    println("Test ScanningStrategy ==> ", ss)
end

@testset "get_pointings-Test" begin    
    pix_tod, psi_tod = get_pointings(ss, 0, 30)
    @test typeof(pix_tod) <: Array
    @test typeof(psi_tod) <: Array
    @show pix_tod[1:10]
    @show psi_tod[1:10]
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

