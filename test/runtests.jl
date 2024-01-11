using Test
include("../src/ME4OH.jl")

@testset "Profile reading" begin
    lon, lat, dates, depth, T, S = read_profile(datafilelist[1])
    @test length(lon) == 8781 
    @test length(lat) == 8781 
    @test lon[1] == 165.55f0
    @test lat[2] == -74.75f0
    @test dates[10] == DateTime(1979, 1, 30)
    @test size(T) == (51, 8781)
    @test T[21] == -1.5671254f0
    @test S[end-100] == 34.910427f0
end