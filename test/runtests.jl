using Falcons
using Healpix
using Test

day = 60 * 60 * 24
year = day * 365
t = 100

function test_wo_imo()
    @testset "get_pointings-Test" begin
        ss = gen_ScanningStrategy()
        pointings = get_pointings(ss, 0, t)
        @test typeof(pointings[1]) <: Matrix
        @test typeof(pointings[2]) <: Matrix
        @test typeof(pointings[3]) <: Matrix
    end
end

function test_w_imo(imo_path)
    ss_imo = gen_ScanningStrategy(duration=t)
    @testset "gen_imo-Test" begin
        imo = gen_imo(imo_path)
        @test typeof(imo)<:Imo
    end

    @testset "get_pointings (IMo ver.) -Test" begin
        imo = gen_imo(imo_path)
        pointings = get_pointings(ss_imo, 0, t)
        @test typeof(pointings[1]) <: Matrix
        @test typeof(pointings[2]) <: Matrix
        @test typeof(pointings[3]) <: Matrix
    end

    @testset "Pick up data from IMo" begin
        imo = gen_imo(imo_path)
        @test_nowarn imo_telescope!(ss_imo, imo, telescope="L")
        @test_nowarn imo_telescope!(ss_imo, imo, telescope="M")
        @test_nowarn imo_telescope!(ss_imo, imo, telescope="H")
        channel_info = get_channel_info(imo)
        for i in channel_info
            @test_nowarn imo_channel!(ss_imo, imo, channel=i)
        end
        bolonames = ["000_000_006_UA_040_T",
            "000_003_008_QB_040_T",
            "000_004_000_UA_040_T",
            "000_007_002_QB_040_T"]
        @test_nowarn imo_name!(ss_imo, imo, name=bolonames)
    end
end

imo_test = Base.prompt("Do you have the schema.json? [y/n]")
if imo_test == "y"
    imo_path = Base.prompt("Input the path for schema.json")
     test_wo_imo()
    test_w_imo(imo_path)
else
    test_wo_imo()
end

test_wo_imo()
