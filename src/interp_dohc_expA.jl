# --------------------------------------------------
# Interpolate the ocean_heat_content_density_anomaly
# --------------------------------------------------
#
# In this script the variables dohc, adohc and dadohc are
# interpolated for the periods of interests. 
#
# A relative error field is also computed for each analysis.

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


# Set variable and period of interest
varname = "adohc"
casename = "fixed-L"
experimentname = "A"
optimize_L = false
optimize_eps = false

thetimeperiodlist = [timeperiod1, timeperiod2, timeperiod3]

# Loop on the 3 time periods
for thetimeperiod in thetimeperiodlist

    @info("Worlking on the time period $(thetimeperiod[1]) - $(thetimeperiod[end])")

    datadirdisk = "$(databasedir)/$(thetimeperiod[1])-$(thetimeperiod[end])"
    isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");
    outputdir = joinpath(mainoutputdir, "experiment-$(experimentname)", varname)
    mkpath(outputdir)

    # Generate list of files for the period of interest 
    datafilelist = ME4OH.get_filelist(datadirdisk);
    timeperiodtext = split(datadirdisk, "/")[end]
    nfiles = length(datafilelist)

    # Create metrics
    _, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)

    @info("Interpolating variable $(varname) over time period $(thetimeperiod[1]) - $(thetimeperiod[end])")
    @info("Experiment $(experimentname)")

    # Loop on the 3 layers
    for (iilayer, layer) in enumerate(depthlayers)
        @info("Working on layer $(layer[1]) - $(layer[end]) m")

        # Create mask for the considered depth levels
        _, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, layer[1])

        # Create empty files for the results
        # (one file per depth)
        fname = ME4OH.make_fname(thetimeperiod, layer, experimentname)
        outputfile = joinpath(outputdir, fname)

        isfile(outputfile) ? rm(outputfile) : @debug("ok") 
        ME4OH.create_netcdf_results(outputfile, varname, longrid, latgrid, thetimeperiod)

        # Loop on the files
        for (iii, datafile) in enumerate(datafilelist)

            @info("Working on file $(basename(datafile)) ($(iii)/$(nfiles))")
            lon, lat, dates, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = ME4OH.read_profile(datafile);

            if varname == "adohc"
                var2interp = adohc
            elseif varname == "dohc"
                var2interp = dohc
            elseif varname == "dadohc"
                var2interp = dadohc
            else 
                @error("Variable $(varname) does not exist")
            end 
            
            goodvalues = .!(isnan.(var2interp[iilayer,:]))

            if optimize_L
                lenxy, infoxy = DIVAnd.fithorzlen((lon[goodvalues], lat[goodvalues], layer[1] * ones(length(lon[goodvalues]))), 
                Float64.(var2interp[1,goodvalues]), layer[1])
                L = (lenxy[1], lenxy[1])
            else
                L = (5., 5.)
            end

            epsilon2 = 10.

            if optimize_eps
                bestfactor_L, bestfactor_eps, cvval, cvvalues, x2Ddata, y2Ddata, cvinter, xi2D, yi2D = 
                DIVAnd_cv(mask, (pm, pn), (xi, yi), (Float64.(lon[goodvalues]), Float64.(lat[goodvalues])),
                Float64.(var2interp[iilayer, goodvalues]), L, epsilon2, 3, 3, 0; moddim=[1,0]);
                L_opt = L .* bestfactor_L
                eps_opt = epsilon2 * bestfactor_eps
                @info("Lopt = $(L_opt), eps_opt = $(eps_opt)")
            end
            
            fi, s = DIVAndrun(mask, (pm, pn), (xi, yi), (Float64.(lon[goodvalues]), Float64.(lat[goodvalues])), 
                Float64.(var2interp[iilayer, goodvalues]), L, epsilon2, moddim=[1,0])
            cpme = DIVAnd_cpme(mask, (pm, pn), (xi, yi), (lon[goodvalues], lat[goodvalues]), 
                var2interp[iilayer, goodvalues], L, epsilon2, moddim=[1,0]);

            NCDataset(outputfile, "a") do ds
                ds[varname][:,:,iii] = fi
                ds[varname*"_error"][:,:,iii] = cpme
            end
            
        end

    end

end
