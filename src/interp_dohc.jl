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

datafilelist = ME4OH.get_filelist(datadirdisk);
timeperiodtext = split(datadirdisk, "/")[end]

# Create empty files for the results
# (one file per depth)

# Create metrics
_, (pm, pn), (xi, yi) = DIVAnd.DIVAnd_rectdom(longrid, latgrid)
_, _, mask = DIVAnd.load_mask(gebco16file, true, longrid, latgrid, 0.)

len_time = zeros(Float64, length(datafilelist))

# Loop on the files
for (iii, datafile) in enumerate(datafilelist)

    @info("Working on file $datafile")
    @time lon, lat, dates, vertical_levels, T, S, dohc, adohc, dadohc, dohc_mask, ts_bounds, depth_level_thickness = 
    ME4OH.read_profile(datafile);

    @info("Fitting the correlation length")
    goodvalues = .!(isnan.(dohc[2,:]))
    lenxy, infoxy = DIVAnd.fithorzlen((lon[goodvalues], lat[goodvalues], zeros(length(lon[goodvalues]))), Float64.(dohc[2,goodvalues]), [0.,])
    len_time[iii] = lenxy[1];

    # DIVAndrun(mask, (pm, pn), (xi, yi), (lon, lat), adohc, lenxy, 5.)
end

# Write the results into a text file
open("correlationlength_$(timeperiodtext)_dohc_level2.txt", "w") do io
    writedlm(io, len_time)
end

