# Compute the Density of Ocean Heat Content from the temperature field
#
# The 4D temperature field is computed by `interp_dohc_expC.jl`, on the regular grid 
# specified by `longrid` and `latgrid` defined in `config.jl`.

using NCDatasets
include("ME4OH.jl")
include("config.jl")    

doplot = false

outputdir = joinpath(mainoutputdir, "experiment-C")
mkpath(outputdir);

thetimeperiod = timeperiod2
layer = depthlayer1

temperaturefile = joinpath(outputdir, "temperature_$(thetimeperiod[1])_$(thetimeperiod[end]).nc")
varname = "dohc"
casename = "fixed-L"
experimentname = "C"

rho0 = 1030    # kg/m^3
Cp0 = 3989.244 # J/(kg K)
T0 = 273.15    # Celsius -> Kelvin
tera = 10^12   # 
scale = Cp0 * rho0/tera

thetimeperiodlist = [timeperiod1, timeperiod2, timeperiod3]

# Loop on the 3 time periods
for thetimeperiod in thetimeperiodlist

    @info("Working on the time period $(thetimeperiod[1]) - $(thetimeperiod[end])")

    # Loop on the 3 layers
    for (iilayer, layer) in enumerate(depthlayers)
        @info("Working on layer $(layer[1]) - $(layer[end]) m")

        # Create a netCDF to store the results
        # (one per layer)
        fname = ME4OH.make_fname(thetimeperiod, layer, "$(experimentname)")
        outputfileprod = joinpath(outputdir, fname)

        # If the file is already there, let's keep it
        isfile(outputfileprod) ? @debug("ok") : ME4OH.create_netcdf_results(outputfileprod, varname, longrid, latgrid, thetimeperiod)
        @info("Writing DOHC in file $(outputfileprod)");

        # Read the climatology and loop over the time periods
        # (loading the full matrix in one go doesn't work)
        NCDataset(temperaturefile, "r") do nc
            depth = nc["depth"][:]
            time = nc["time"][:]
            gooddepth = findall((depth .>= depthlayer1[1]) .& (depth .<= depthlayer1[2]))
            # push!(gooddepth, gooddepth[end] + 1)
            #@info(length(gooddepth))

            Δdepth = diff(depth)[gooddepth]
            Δdepth3D = repeat(Δdepth, outer=(1, length(longrid), length(latgrid)))
            Δdepth3D = permutedims(Δdepth3D, [2,3,1])
            #@show(size(Δdepth3D))

            # Loop on time
            for (ii, tt) in enumerate(time)
                T = nc["temperature"][:,:,gooddepth, ii]
                ok = (T .+ T0) .* Δdepth3D;
                dohc = scale .* sum(ok, dims=3);
                dohc = dropdims(dohc; dims=3)
                #@show(size(dohc));
                
                if doplot
                    fig = plt.figure(figsize=(12, 8))
                    ax = plt.subplot(111)
                    pcm = ax.pcolormesh(longrid, latgrid, dohc')
                    cb = plt.colorbar(pcm)
                    plt.savefig(joinpath(figdir, "test_dohc.jpg"))
                    plt.close()
                end

                NCDataset(outputfileprod, "a") do ds
                    ds[varname][:,:,ii] = dohc
                end

            end

        end
    end
end