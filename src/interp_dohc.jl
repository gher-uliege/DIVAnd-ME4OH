# Interpolate the ocean_heat_content_density_anomaly
# --------------------------------------------------

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

experimentname = "A"
thetimeperiod = timeperiod1
datadirdisk = "/media/ctroupin/T7 Shield/000060826v009/data/en4.1.1/$(thetimeperiod[1])-$(thetimeperiod[end])"
isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");

datafilelist = ME4OH.get_filelist(datadirdisk);
timeperiodtext = split(datadirdisk, "/")[end]

# Create metrics
_, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)
_, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, 0.)

# Loop on the 3 layers
for (iilayer, layer) in enumerate(depthlayers)
    @info("Working on layer $(layer[1]) - $(layer[end]) m")

    # Create empty files for the results
    # (one file per depth)
    fname = ME4OH.make_fname(thetimeperiod, layer, experimentname)
    outputfile = joinpath(outputdir, fname)


    isfile(outputfile) ? rm(outputfile) : @debug("ok") 
    ME4OH.create_netcdf_results(outputfile, longrid, latgrid, thetimeperiod)

    # Loop on the files
    for (iii, datafile) in enumerate(datafilelist)

        @info("Working on file $datafile")
        @time lon, lat, dates, vertical_levels, T, S, dohc, adohc, dadohc, dohc_mask, ts_bounds, depth_level_thickness = 
        ME4OH.read_profile(datafile);
        goodvalues = .!(isnan.(adohc[iilayer,:]))

        lenxy, infoxy = DIVAnd.fithorzlen((lon[goodvalues], lat[goodvalues], zeros(length(lon[goodvalues]))), 
        Float64.(adohc[1,goodvalues]), layer[1])
    
        fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (Float64.(lon[goodvalues]), Float64.(lat[goodvalues])), 
            Float64.(adohc[iilayer, goodvalues]), (lenxy[1], lenxy[1]), 10., moddim=[1,0])

        NCDataset(outputfile, "a") do ds
            ds["dohc"][:,:,iii] = fi
        end

    end

end
