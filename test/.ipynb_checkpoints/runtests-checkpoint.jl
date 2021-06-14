using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
prec1 = 60 * 60 * 3


nside = 128
times = year #[sec]
sampling_rate = 1 #[Hz]
FP_theta = [0.0] #The angle with respect to the boresight, 0 degree represents the boresight.
FP_phi = [0.0]
alpha = 55.0 #[degree]
beta = 60.0 #[degree]
prec_period = 180.22 #[min]
spin_rpm = 0.04 #[rpm]
hwp_rpm = 0.05 #[rpm]
start_point = "pole" #You can choose "pole" or "equator"

ss = ScanningStrategy(
    nside,
    times,
    sampling_rate,
    alpha,
    beta,
    prec_period,
    spin_rpm,
    hwp_rpm,
    FP_theta,
    FP_phi,
    start_point)

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

