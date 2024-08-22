# Construct a monthly climatology for Temperature
# -----------------------------------------------

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
experimentname = "C-SST"
optimize_L = false
optimize_eps = false

# analysis parameters 
Lbackground = (15., 15.)
Lproduct = (5., 5.)
eps2background = 10.
eps2product = 5.

thetimeperiod = timeperiod2
datadirdisk = joinpath(databasedir, "$(thetimeperiod[1])-$(thetimeperiod[end])")
isdir(datadirdisk) ? @info("Directory exists") : @error("Directory does not exist");
outputdir = joinpath(mainoutputdir, "experiment-C")
mkpath(outputdir)

timeperiodtext = "$(thetimeperiod[1])01_$(thetimeperiod[end])12"

fnameprefix = "ofam3-jra55.all.EN.4.1.1.f.profiles.g10"
fnamesuffix = "update.extra.danom"

# Create metrics
_, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)

outputdir = joinpath(mainoutputdir, "experiment-$(experimentname)", casename)
mkpath(outputdir)

ntimes = length(thetimeperiod) * 12

@info("Creating gridded field for the temperature")

# Get depth levels and level thickness from the 1st file
datafilelist = ME4OH.get_filelist(datadirdisk);
_, _ , _, depthlevels, _, _, depth_level_thickness = ME4OH.read_TS_data(datafilelist[1])
depthlevels = unique(depthlevels)

# Create file for the climatology
fname = "temperature_climatology_$(thetimeperiod[1])_$(thetimeperiod[end]).nc"
outputfileclim = joinpath(outputdir, fname)
isfile(outputfileclim) ? rm(outputfileclim) : @debug("ok") 
ME4OH.create_netcdf_climatology4D(outputfileclim, longrid, latgrid, depthlevels, thetimeperiod)

# Create file for the product
fname = "temperature_$(thetimeperiod[1])_$(thetimeperiod[end]).nc"
outputfileprod = joinpath(outputdir, fname)
isfile(outputfileprod) ? rm(outputfileprod) : @debug("ok") 
ME4OH.create_netcdf_product4D(outputfileprod, longrid, latgrid, depthlevels, thetimeperiod)

# Loop on the months
for mm = 1:12

    @info("Working on month $(mm) / 12")
    
    # Find all the files for the considered month
    datafilelistmonth = glob("$(fnameprefix).*$(lpad(mm, 2, "0")).$(fnamesuffix).$(timeperiodtext).nc", datadirdisk)
    @info("Found $(length(datafilelistmonth)) files")

    # Load all the data for that month
    obslon, obslat, obsdates, obsdepth, T, S = ME4OH.read_TS_data(datafilelistmonth);
    nobs = length(obslon)
    @info("Total number of observations: $(nobs)")
    
    # Loop on the depths
    for (iilayer, layer) in enumerate(depthlevels)
        @info("Working on depth: $(layer) m")

        # Create mask for the considered depth levels
        _, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, layer)

        # Select the data for the corresponding depth
        gooddepth = findall(obsdepth .== layer)
        @info("Found $(length(gooddepth)) observations")

        obslonlayer = @view(obslon[gooddepth])
        obslatlayer = @view(obslat[gooddepth])
        obsdateslayer = @view(obsdates[gooddepth])
        Tlayer = @view(T[gooddepth])

        # Perform analysis with all the coordinates for that month
        ficlim, s = DIVAndrun(mask, (pm, pn), (xi, yi), (obslon[gooddepth], obslat[gooddepth]), 
            T[gooddepth], Lbackground, eps2background, moddim=[1,0])

        # Compute associated error field
        cpmeclim = DIVAnd_cpme(mask, (pm, pn), (xi, yi), (obslonlayer, obslatlayer), 
            Tlayer, Lbackground, eps2background, moddim=[1,0]);

        if doplot
            fig = plt.figure(figsize=(10, 8))
            ax = plt.subplot(111)
            ax.pcolormesh(longrid, latgrid, ficlim', cmap=plt.cm.RdYlBu_r)
            fname = "temperature_clim_$(lpad(mm, 2, '0'))_depth_$(iilayer).jpg"
            ax.set_title("Temperature at $(layer) m for month $(Dates.monthname(mm))")
            plt.savefig(joinpath(figdir, fname))
            plt.close()
        end

        # Write the results in the climatology file
        NCDataset(outputfileclim, "a") do ds
            ds["temperature"][:,:,iilayer,mm] = ficlim
            ds["temperature_error"][:,:,iilayer,mm] = cpmeclim
        end

        # Compute the residuals (from the results)
        # (they correspond to a given layer and given month)
        dataresidual = DIVAnd_residual(s, ficlim);

        # Start loop on the years 
        # (according to the selected period)

        for yyyy in thetimeperiod
            @info("Working on years: $(yyyy)")

            # Selec the observations for that year
            yearsel = findall(Dates.year.(obsdateslayer) .== yyyy)
            @info("Found $(length(yearsel)) observations")

            # Perform analysis on the anomalies
            fimonth, s = DIVAndrun(mask, (pm, pn), (xi, yi), (obslonlayer[yearsel], obslatlayer[yearsel]), 
                dataresidual[yearsel], Lproduct, eps2product, moddim=[1,0])

            # Compute associated error field
            cpme = DIVAnd_cpme(mask, (pm, pn), (xi, yi), (obslonlayer[yearsel], obslatlayer[yearsel]), 
                dataresidual[yearsel], Lproduct, eps2product, moddim=[1,0]);

            field2plot = fimonth .+ ficlim

            if doplot
                fig = plt.figure(figsize=(10, 8))
                ax = plt.subplot(111)
                ax.pcolormesh(longrid, latgrid, field2plot', cmap=plt.cm.RdYlBu_r)
                fname = "temperature_$(yyyy)_$(lpad(mm, 2, '0'))_depth_$(iilayer).jpg"
                ax.set_title("Temperature at $(layer) m for month $(Dates.monthname(mm)) and year $(yyyy)")
                plt.savefig(joinpath(figdir, fname))
                plt.close()
            end

            # Write the results in the netCDF file
            timeindex = (yyyy - thetimeperiod[1]) * 12 + mm  

            @info("Writing time step $(timeindex) / $(ntimes)")
            NCDataset(outputfileprod, "a") do ds
                ds["temperature"][:,:,iilayer,timeindex] = field2plot
                ds["temperature_error"][:,:,iilayer,timeindex] = cpme
            end

        end

    end

end

    # # Create empty files for the results
    # # (one file per depth)
    # fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)")
    # outputfileprod = joinpath(outputdir, fname)
    # #isfile(outputfileprod) ? rm(outputfileprod) : @debug("ok") 
    # #ME4OH.create_netcdf_results(outputfileprod, "dohc", longrid, latgrid, thetimeperiod)

    # # Create empty files for the anomalies
    # # (to compare with experiment A)
    # fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)-anom")
    # outputfileanom = joinpath(outputdir, fname)
    # isfile(outputfileanom) ? rm(outputfileanom) : @debug("ok") 
    # ME4OH.create_netcdf_results(outputfileanom, "adohc", longrid, latgrid, thetimeperiod)

    

#         
        
#         if doplot
#             fig = plt.figure(figsize=(10, 8))
#             plt.pcolormesh(longrid, latgrid, ficlim', cmap=plt.cm.RdYlBu_r, vmin=0.34, vmax=0.38)
#             fname = "dohc_clim_$(lpad(mm, 2, '0'))_depth_$(iilayer).jpg"
#             plt.savefig(joinpath(figdir, fname))
#             plt.close()
#         end

#         
#        

#     end

# end
