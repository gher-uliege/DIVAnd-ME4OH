# Construct a monthly climatology for DOHC
# ----------------------------------------

using DIVAnd
using Dates
using NCDatasets
using Test
using PyPlot
using Statistics
using DelimitedFiles
using Glob
const plt = PyPlot
include("./ME4OH.jl")
include("./config.jl")

doplot = false
varname = "dohc"
casename = "fixed-L"
experimentname = "B-mask"
optimize_L = false
optimize_eps = false
maskdepth = true

# analysis parameters 
Lbackground = (15., 15.)
Lproduct = (5., 5.)
eps2background = 10.
eps2product = 5.

thetimeperiod = timeperiod2
datadirdisk = joinpath(databasedir, "$(thetimeperiod[1])-$(thetimeperiod[end])")
isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");
outputdir = joinpath(mainoutputdir, "experiment-B")
mkpath(outputdir)

# timeperiodtext = split(datadirdisk, "/")[end]
timeperiodtext = "$(thetimeperiod[1])01_$(thetimeperiod[end])12"

fnameprefix = "ofam3-jra55.all.EN.4.1.1.f.profiles.g10"
fnamesuffix = "update.extra.danom"

# Create metrics
_, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)

outputdir = joinpath(mainoutputdir, "experiment-$(experimentname)", casename)
mkpath(outputdir)

ntimes = length(thetimeperiod) * 12

@info("Creating gridded field for the Ocean Heat Content Density")

# Loop on the 3 layers
for (iilayer, layer) in enumerate(depthlayers)
    @info("Working on layer $(layer[1]) - $(layer[end]) m")

    # Create mask for the considered depth levels
    _, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, layer[1])

    # Create empty files for the results
    # (one file per depth)
    fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)")
    outputfileprod = joinpath(outputdir, fname)
    isfile(outputfileprod) ? rm(outputfileprod) : @debug("ok") 
    ME4OH.create_netcdf_results(outputfileprod, "dohc", longrid, latgrid, thetimeperiod)

    # Create a netCDF for the climatology
    fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)-clim")
    outputfileclim = joinpath(outputdir, fname)
    isfile(outputfileclim) ? rm(outputfileclim) : @debug("ok") 
    ME4OH.create_netcdf_climatology(outputfileclim, "dohc", longrid, latgrid, thetimeperiod)

    # Create empty files for the anomalies
    # (to compare with experiment A)
    fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)-anom")
    outputfileanom = joinpath(outputdir, fname)
    isfile(outputfileanom) ? rm(outputfileanom) : @debug("ok") 
    ME4OH.create_netcdf_results(outputfileanom, "adohc", longrid, latgrid, thetimeperiod)

    # Loop on the months
    for mm = 1:12

        @info("Working on month $(mm) / 12")
        
        # Find all the files for the considered month
        datafilelistmonth = glob("$(fnameprefix).*$(lpad(mm, 2, "0")).$(fnamesuffix).$(timeperiodtext).nc", datadirdisk)
        @info("Found $(length(datafilelistmonth)) files")

        # Load all the data for that month
        lon, lat, dates, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = ME4OH.read_profile(datafilelistmonth);
        nobs = length(lon)
        @info("Total number of observations: $(nobs)")
        
        # Remove NaNs
        goodvalues = .!(isnan.(dohc[iilayer,:]))
        ngood = sum(goodvalues)
        pct = round(100 * ngood/nobs, digits=2)
        @info("Number of good observations: $(sum(goodvalues))/$(nobs) ($(pct)%)")

        if maskdepth == true
            # Use the variable "dohc_mask_by_en4_maxdepth" to select valid obs.
            goodmask = dohc_mask[iilayer,:] .== 1
            @info("Number of non-masked observations: $(sum(goodmask))")
            goodvalues = goodvalues .& goodmask
            @info("Final number of good observations: $(sum(goodvalues))")
        end

        obslon = Float64.(lon[goodvalues]) 
        obslat = Float64.(lat[goodvalues])
        obsdates = dates[goodvalues]
        obsval = Float64.(dohc[iilayer, goodvalues])

        # Perform analysis with all the coordinates for that month
        ficlim, s = DIVAndrun(mask, (pm, pn), (xi, yi),(obslon, obslat), obsval, 
            Lbackground, eps2background, moddim=[1,0])
        # Compute associated error field
        cpme = DIVAnd_cpme(mask, (pm, pn), (xi, yi),(obslon, obslat), obsval, 
            Lbackground, eps2background, moddim=[1,0]);
        
        if doplot
            fig = plt.figure(figsize=(10, 8))
            plt.pcolormesh(longrid, latgrid, ficlim', cmap=plt.cm.RdYlBu_r, vmin=0.34, vmax=0.38)
            fname = "dohc_clim_$(lpad(mm, 2, '0'))_depth_$(iilayer).jpg"
            plt.savefig(joinpath(figdir, fname))
            plt.close()
        end

        # Write the results in the climatology file
        NCDataset(outputfileclim, "a") do ds
            ds["dohc"][:,:,mm] = ficlim
            ds["dohc_error"][:,:,mm] = cpme
        end

        # Compute the residuals (from the results)
        dataresidual = DIVAnd_residual(s, ficlim);

        # Start loop on the years 
        # (according to the selected period)
        for yyyy in thetimeperiod
            @info("Working on years: $(yyyy)")

            # Selec the observations for that year
            yearsel = findall(Dates.year.(obsdates) .== yyyy)
            @info("Found $(length(yearsel)) observations")

            # Perform analysis on the anomalies
            @time fimonth, s = DIVAndrun(mask, (pm, pn), (xi, yi), (obslon[yearsel], obslat[yearsel]), 
            dataresidual[yearsel], Lproduct, eps2product, moddim=[1,0])

            # Compute associated error field
            cpme = DIVAnd_cpme(mask, (pm, pn), (xi, yi), (obslon[yearsel], obslat[yearsel]), 
            dataresidual[yearsel], Lproduct, eps2product, moddim=[1,0]);

            field2plot = fimonth .+ ficlim

            if doplot
                fig = plt.figure(figsize=(10, 8))
                plt.pcolormesh(longrid, latgrid, field2plot', cmap=plt.cm.RdYlBu_r, vmin=0.34, vmax=0.38)
                fname = "dohc_$(yyyy)_$(lpad(mm, 2, '0'))_depth_$(iilayer).jpg"
                plt.savefig(joinpath(figdir, fname))
                plt.close()
            end

            # Write the results in the netCDF file
            timeindex = (yyyy - thetimeperiod[1]) * 12 + mm  

            @info("Writing time step $(timeindex) / $(ntimes)")
            NCDataset(outputfileprod, "a") do ds
                ds["dohc"][:,:,timeindex] = field2plot
                ds["dohc_error"][:,:,timeindex] = cpme
            end

            NCDataset(outputfileanom, "a") do ds
                ds["adohc"][:,:,timeindex] = fimonth
                ds["adohc_error"][:,:,timeindex] = cpme
            end

        end

    end

end
