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
experimentname = "A-time"
optimize_L = false
optimize_eps = false
Ltime = 30.                     # temporal correlation length, in month(s) 
nmonth = 4                      # number of months to consider before and after
dateref = Dates.DateTime(1970,1,1)

thetimeperiodlist = [timeperiod1, timeperiod2, timeperiod3]

# Loop on the 3 time periods
for thetimeperiod in thetimeperiodlist

    @info("Working on the time period $(thetimeperiod[1]) - $(thetimeperiod[end])")

    datadirdisk = joinpath(databasedir, "$(thetimeperiod[1])-$(thetimeperiod[end])")
    isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");
    outputdir = joinpath(mainoutputdir, "experiment-$(experimentname)", varname)
    mkpath(outputdir)

    # Generate list of files for the period of interest 
    datafilelist = ME4OH.get_filelist(datadirdisk);
    timeperiodtext = split(datadirdisk, "/")[end]
    nfiles = length(datafilelist)

    @info("Interpolating variable $(varname) over time period $(thetimeperiod[1]) - $(thetimeperiod[end])")
    @info("Experiment $(experimentname)")

    # Loop on the 3 layers
    for (iilayer, layer) in enumerate(depthlayers)
        @info("Working on layer $(layer[1]) - $(layer[end]) m")

        # Create mask for the considered depth levels
        _, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, layer[1])

        # Create empty files for the results
        # (one file per depth)
        fname = ME4OH.make_fname(thetimeperiod, layer, "A")
        outputfile = joinpath(outputdir, fname)

        isfile(outputfile) ? rm(outputfile) : @debug("ok") 
        ME4OH.create_netcdf_results(outputfile, varname, longrid, latgrid, thetimeperiod, title="DIVAnd interpolated field with temporal correlation")

        # Loop on the files
        for (iii, datafile) in enumerate(datafilelist)

            # Load central time
            lon, lat, dates, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = ME4OH.read_profile(datafile);
            datecentral = Dates.DateTime(Dates.year(dates[1]), Dates.month(dates[1]), 15)
            dategrid = datecentral - Dates.Month(nmonth) : Dates.Month(1) : datecentral + Dates.Month(nmonth)
            timegrid = Dates.value.(dategrid - dateref) / (86400000)
            tindex = findfirst(dategrid.==datecentral)

            _, (pm, pn, pt), (xi, yi, ti) = DIVAnd.DIVAnd_rectdom(longrid, latgrid, timegrid)

            # Repeat the mask along time dimension
            masktime = repeat(mask, outer=(1,1,length(timegrid)))

            tmin = maximum([1, iii - nmonth])
            tmax = minimum([nfiles, iii + nmonth])

            @info("Working on files before & after $(basename(datafile)) [$(tmin)-$(tmax)] ($(iii)/$(nfiles))")
            lon, lat, dates, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = ME4OH.read_profile(datafilelist[tmin:tmax]);

            obsdays = Dates.value.(dates - dateref) / (86400 * 1000)
            @info(extrema(obsdays))
            @info(extrema(timegrid))

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
            @info("Number of good observations: $(sum(goodvalues))");
            @info(size(var2interp))

            if optimize_L
                lenxy, infoxy = DIVAnd.fithorzlen((lon[goodvalues], lat[goodvalues], layer[1] * ones(length(lon[goodvalues]))), 
                Float64.(var2interp[1,goodvalues]), layer[1])
                L = (lenxy[1], lenxy[1])
            else
                L = (5., 5., 30.)
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
            
            @info("Creating gridded and error fields")
            @time fi, s = DIVAndrun(masktime, (pm, pn, pt), (xi, yi, ti), (lon[goodvalues], lat[goodvalues], obsdays[goodvalues]), 
                var2interp[iilayer, goodvalues], L, epsilon2, moddim=[1,0,0])
            cpme = DIVAnd_cpme(masktime, (pm, pn, pt), (xi, yi, ti), (lon[goodvalues], lat[goodvalues], obsdays[goodvalues]), 
                var2interp[iilayer, goodvalues], L, epsilon2, moddim=[1,0,0]);

            @info(size(fi));

            NCDataset(outputfile, "a") do ds
                ds[varname][:,:,iii] = fi[:,:,tindex]
                ds[varname*"_error"][:,:,iii] = cpme[:,:,tindex]
            end
            
        end

    end

end
