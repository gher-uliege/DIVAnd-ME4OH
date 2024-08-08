# Construct a monthly climatology for DOHC
# ----------------------------------------

using DIVAnd
using Dates
using NCDatasets
using Test
using PyPlot
using Statistics
using DelimitedFiles
const plt = PyPlot
include("./ME4OH.jl")
include("./config.jl")

varname = "dohc"
casename = "fixed-L"
experimentname = "B"
optimize_L = false
optimize_eps = false

L = (15., 15.)
epsilon2 = 10.

thetimeperiod = timeperiod1
datadirdisk = "/media/ctroupin/T7 Shield/000060826v009/data/en4.1.1/$(thetimeperiod[1])-$(thetimeperiod[end])"
isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");
timeperiodtext = split(datadirdisk, "/")[end]

# Create metrics
_, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)

outputdir = joinpath(mainoutputdir, "experiment-$(experimentname)", "climatology", casename)
mkpath(outputdir)

# Loop on the 3 layers
for (iilayer, layer) in enumerate(depthlayers)
    @info("Working on layer $(layer[1]) - $(layer[end]) m")

    # Create mask for the considered depth levels
    _, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, layer[1])

    # Create empty files for the results
    # (one file per depth)
    fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)-clim")
    outputfile = joinpath(outputdir, fname)

    isfile(outputfile) ? rm(outputfile) : @debug("ok") 
    #ME4OH.create_netcdf_results(outputfile, varname, longrid, latgrid, thetimeperiod)

    # Loop on the month
    for mm = 1:12

        @info("Working on month $(mm)/12)")
        

        datafilelistmonth = glob("ofam3-jra55.all.EN.4.1.1.f.profiles.g10.*$(lpad(mm, 2, "0")).update.extra.danom.197901_201412.nc", datadirdisk)
        @info("Found $(length(datafilelistmonth)) files")

        lon, lat, dates, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = ME4OH.read_profile(datafilelistmonth);

        @info("Number of observations: $(length(lon))")
        
        goodvalues = .!(isnan.(dohc[iilayer,:]))

        #fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (Float64.(lon[goodvalues]), Float64.(lat[goodvalues])), 
        #    Float64.(dohc[iilayer, goodvalues]), L, epsilon2, moddim=[1,0])

        #@time lon, lat, dates, vertical_levels, T, S, dohc, adohc, dadohc, dohc_mask, ts_bounds, depth_level_thickness = 
        #ME4OH.read_profile(datafile);

    end

end
