using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
prec1 = 60 * 60 * 3

ss = scan_strategy()
ss.nside = 128
ss.times = day
ss.sampling_rate = 2
ss.FP_theta = [0.0]
ss.FP_phi = [0.0]
ss.alpha = 55
ss.beta = 60
ss.prec_period = 180.223
ss.spin_rpm = 0.04
ss.hwp_rpm = 0.05
ss.start_point = "pole"

@testset "scan_strategy_structure-Test" begin
    @show ss.nside
    @show ss.times
    @show ss.sampling_rate
    @show ss.FP_theta
    @show ss.FP_phi
    @show ss.alpha
    @show ss.beta
    @show ss.prec_period
    @show ss.spin_rpm
    @show ss.hwp_rpm
    @test typeof(ss) <: scan_strategy
end

@testset "get_scan_tod-Test" begin    
    pix_tod, psi_tod = get_scan_tod(ss, 0, 30)
    @show pix_tod[1:10]
    @show psi_tod[1:10]
    @test typeof(pix_tod) <: Array 
    @test typeof(psi_tod) <: Array 
end

@testset "Mapmaking-Test" begin
    n = 2
    outmap = Mapmaking(ss, 2)
    @test length(outmap[1]) == nside2npix(ss.nside)
    @show outmap[1][1:5]
    @show outmap[2][1:5]
    @show outmap[3][1:5]
    @show outmap[4][1:5]
    @show outmap[5][1:5]
end

