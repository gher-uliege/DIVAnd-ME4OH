using Test
using Dates

include("../src/config.jl")
include("../src/ME4OH.jl")

@testset "listing files" begin
    download_check(datafiletest1, datafiletest1URL)
    download_check(datafiletest2, datafiletest2URL)
    download_check(datafiletest3, datafiletest3URL)

    datafilelist1 = ME4OH.get_filelist(datatestdir, 1970:2000)
    datafilelist2 = ME4OH.get_filelist(datatestdir, 1990:2000)
    @test length(datafilelist1) == 3
    @test length(datafilelist2) == 0

end

@testset "reading profile" begin
    download_check(datafiletest1, datafiletest1URL)
    download_check(datafiletest3, datafiletest3URL)

    lon, lat, dates, depth, T, S, dohc, dohc_mask = ME4OH.read_profile(datafiletest1)
    @test length(lon) == 8781 
    @test length(lat) == 8781 
    @test lon[1] == 165.55f0
    @test lat[2] == -74.75f0
    @test dates[10] == DateTime(1979, 1, 30)
    @test size(T) == (51, 8781)
    @test T[21] == -1.5671254f0
    @test S[end-100] == 34.910427f0
    @test dohc[3, 10] == 1.5131863f0
    @test sum(isnan.(dohc)) == 7708
    @test sum(dohc_mask) == length(lon) * 3

    lon, lat, dates, depth, T, S, dohc, dohc_mask = ME4OH.read_profile(datafiletest3);
    @test length(lon) == 8781 
    @test length(lat) == 8781 
    @test lon[1] == 165.55f0
    @test lat[2] == -74.75f0
    @test dates[10] == DateTime(1979, 1, 30)
    @test size(T) == (51, 8781)
    @test T[21] == -1.5671254f0
    @test S[end-100] == 34.910427f0
    @test dohc[3, 10] == 1.5131863f0
    @test sum(isnan.(dohc)) == 7708
        @test sum(dohc_mask) == 5398

end

@testset "read data from a file list" begin
    download_check(datafiletest1, datafiletest1URL)
    download_check(datafiletest2, datafiletest2URL)
    datafilelisttest = [datafiletest1, datafiletest2]
    obslon, obslat, obsdepth, obsdates, T, S = ME4OH.read_data(datafilelisttest)

    @test length(obslon) == 767726
    @test length(obslon) == length(obslat)
    @test length(obslon) == length(obsdepth)
    @test length(obslon) == length(obsdates)

    obslon[400] == 192.95f0
    obslat[400] == -71.45f0
    obsdepth[400] == 4509.18f0
    obsdates[400] == DateTime(1979, 01, 30)
    T[end] == 1.6556292f0
    S[end] == 34.979248f0

end

@testset "get dates" begin
    download_check(datafiletest1, datafiletest1URL)
    year = ME4OH.get_year(datafiletest1)
    @test year == 1979
end

@testset "vectorize" begin
    download_check(datafiletest1, datafiletest1URL)

    lon, lat, dates, depth, T, S, dohc = ME4OH.read_profile(datafiletest1)
    obslon, obslat, obsdates, obsdepth, T, S = ME4OH.vectorize_obs(lon, lat, dates, depth, T, S)
    
    @test length(obslon) == 354170
    @test sum(isnan.(obslat)) == 0
    @test obslat[22] == -74.75f0
    @test obsdepth[end] == 195.0f0
    @test obsdates[end-100] == DateTime(1979, 1, 14)
    
end

@testset "create file name" begin
    timeperiod1 = 1979:2014
    depthlayer3 = [685.9, 1985.3] 
    fname1 = ME4OH.make_fname(timeperiod1, depthlayer3, "A")
    @test fname1 == "OHC_1979_2014_lev0_1985.3_expA_DIVAnd.nc"
    fname2 = ME4OH.make_fname(timeperiod1, depthlayer3, "A"; product="Test")
    @test fname2 == "OHC_1979_2014_lev0_1985.3_expA_Test.nc"

end

@testset "getting time period" begin
    timeperiodtest = 2005:2014
    timegrid = ME4OH.get_timegrid(timeperiodtest)
    @test length(timegrid) == 120
    @test timegrid[50] == DateTime(2009, 2, 1)

end

@testset "converting DateTime to days" begin
    timeperiodtest = 2005:2020
    daygrid = ME4OH.datetime2days(timeperiodtest)
    @test length(daygrid) == 192
    @test daygrid[1] == 38351.0

    daygrid2 = ME4OH.datetime2days(timeperiodtest, dateref = DateTime(1980, 1, 1))
    @test length(daygrid2) == 192
    @test daygrid2[1] == 9132.0
end

@testset "counting profiles" begin
    datafilelist = [joinpath(datatestdir, "ofam3-jra55.all.EN.4.1.1.f.profiles.g10.197901.update.nc")]
    nprof = ME4OH.get_profile_number(datafilelist)
    @test length(datafilelist) == length(nprof)
    @test nprof == [8781,]
end