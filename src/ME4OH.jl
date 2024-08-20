module ME4OH

using NCDatasets
using DataStructures
using Dates
using Glob 

"""
    get_filelist(datadir, timeperiod)

Return the list of files lcated in `datadir` that match the time period defined by `timeperiod` (range).

__Notes:__ 
- the files with the `nonan` suffix are not included in the list;
- only the files with the detrended anomalies (`danom`) are included.

# Example
```julia-repl
julia> get_filelist(datadir, 1995:2005)
72-element Vector{String}:
 "/media/ctroupin/T7 Shield/00006" ⋯ 83 bytes ⋯ "te.extra.danom.197901_201412.nc"
 "/media/ctroupin/T7 Shield/00006" ⋯ 83 bytes ⋯ "te.extra.danom.197901_201412.nc"
 ⋮
 "/media/ctroupin/T7 Shield/00006" ⋯ 83 bytes ⋯ "te.extra.danom.197901_201412.nc"
 "/media/ctroupin/T7 Shield/00006" ⋯ 83 bytes ⋯ "te.extra.danom.197901_201412.nc"
```
"""
function get_filelist(datadir::AbstractString, timeperiod::UnitRange{Int64}=1900:2100)::Vector{String}
    datafilelist = Glob.glob("ofam3-jra55.all.EN.4.1.1.f.profiles.g10.*.update.extra.danom.*.nc", datadir)

    # Don't take the nonan files
    filter!(x -> !occursin("nonan", x), datafilelist)
    yearlist = get_year.(datafilelist)
    # Select the good year period
    goodyear = findall((yearlist .<= timeperiod[end]) .& (yearlist .>= timeperiod[1]))
    @info("Found $(length(datafilelist[goodyear])) files")
    return datafilelist[goodyear]
end

"""
    read_profile(datafile)

Read the coordinates and the measurements from a data file

# Example
```julia-repl 
julia> lon, lat, dates, vertical_levels, T, S, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = read_profile(datafile)
```
"""
function read_profile(datafile::AbstractString)
    NCDataset(datafile, "r") do ds
        lon = Float64.(ds["ts_lon"][:])
        lon[lon .< 20.] .+= 360;

        lat = Float64.(ds["ts_lat"][:])
        time = ds["en4_ymd"][:,:]
        dates = Dates.DateTime.(time[1,:], time[2,:], time[3,:])
        vertical_levels = ds["ts_z"][:]
        depth_level_thickness = Float64.(ds["ts_dz"][:])
        bounds = Float64.(ds["ts_bound"][:,:])
        #T = Float64.(ds["temp"][:,:])
        #SSH = Float64.(ds["eta_t"][:])
        #SST = Float64.(ds["sst"][:])
        #S = Float64.(ds["salt"][:,:])
        dohc = Float64.(ds["dohc"][:,:])
        adohc = Float64.(ds["adohc"][:,:])
        dadohc = Float64.(ds["dadohc"][:,:])

        varlist = keys(ds)
        if "dohc_mask_by_en4_maxdepth" in varlist
            @debug("Reading dohc_mask_by_en4_maxdepth from file")
            dohc_mask = Bool.(ds["dohc_mask_by_en4_maxdepth"][:,:])
        else
            dohc_mask = ones(Bool, size(dohc))
        end

        return lon::Vector{Float64}, lat::Vector{Float64}, dates::Vector{DateTime}, 
        # vertical_levels::Vector{Float32}, #T::Matrix{Float32}, S::Matrix{Float32}, 
        dohc::Matrix{Float64}, adohc::Matrix{Float64}, dadohc::Matrix{Float64},
        dohc_mask::AbstractArray{Bool}, bounds::Matrix{Float64}, depth_level_thickness::Vector{Float64}
    end
end

"""
    read_profile(datafilelist)

Read the coordinates and the measurements from a list of files

# Example
```julia-repl 
julia> lon, lat, dates, vertical_levels, T, S, dohc, adohc, dadohc, dohc_mask, bounds, depth_level_thickness = read_profile(datafilelistmonth)
```
"""
function read_profile(datafilelist::Vector{String})

    # Create empty vectors and matrices
    lonall = Float64[]
    latall = Float64[]
    datesall = DateTime[]
    dohcall = Array{Float64}(undef, 3, 0)
    adohcall = Array{Float64}(undef, 3, 0)
    dadohcall = Array{Float64}(undef, 3, 0)
    dohc_maskall = Array{Bool}(undef, 3, 0)

    _, _, _, _, _, _, _, bounds, depth_level_thickness = read_profile(datafilelist[1])

    for datafile in datafilelist

        lon, lat, dates, dohc, adohc, dadohc, dohc_mask, _, _ = read_profile(datafile);

        append!(lonall, lon)
        append!(latall, lat)
        append!(datesall, dates)
        dohcall = hcat(dohcall, dohc)
        adohcall = hcat(adohcall, adohc)
        dadohcall = hcat(dadohcall, dadohc)
        
    end

    return lonall::Vector{Float64}, latall::Vector{Float64}, datesall::Vector{DateTime},
        dohcall::Matrix{Float64}, adohcall::Matrix{Float64}, dadohcall::Matrix{Float64},
        dohc_maskall::AbstractArray{Bool}, bounds::Matrix{Float64}, depth_level_thickness::Vector{Float64}
end

"""
    read_TS_data(datafile)

Read the coordinates and the temperature and salinity measurements from a data file

# Example
```julia-repl 
julia> lon, lat, dates, depth, T, S, depth_level_thickness = read_TS_data(datafile)
```
"""
function read_TS_data(datafile::AbstractString)
    NCDataset(datafile, "r") do ds
        lon = Float64.(ds["ts_lon"][:])
        lat = Float64.(ds["ts_lat"][:])
        time = ds["en4_ymd"][:,:]
        dates = Dates.DateTime.(time[1,:], time[2,:], time[3,:])
        depth = Float64.(ds["ts_z"][:])
        depth_level_thickness = Float64.(ds["ts_dz"][:])
        T = Float64.(ds["temp"][:,:])
        S = Float64.(ds["salt"][:,:])
   
        varlist = keys(ds)
        if "dohc_mask_by_en4_maxdepth" in varlist
            @debug("Reading dohc_mask_by_en4_maxdepth from file")
            dohc_mask = Bool.(ds["dohc_mask_by_en4_maxdepth"][:,:])
        else
            dohc_mask = ones(Bool, size(dohc))
        end

        obslon, obslat, obsdates, obsdepth, obsval1, obsval2 = vectorize_obs(lon, lat, dates, depth, T, S)

        return obslon::Vector{Float64}, obslat::Vector{Float64}, obsdates::Vector{DateTime}, 
        obsdepth::Vector{Float64}, obsval1::Vector{Float64}, obsval2::Vector{Float64}, 
        depth_level_thickness::Vector{Float64}
    end
end

"""
    read_TS_data(datafilelist)

Read the coordinates and the temperature and salinity measurements from a list of data files

# Example
```julia-repl 
julia> lon, lat, dates, depth, T, S, depth_level_thickness = read_TS_data(datafilelist)
```
"""
function read_TS_data(datafilelist::Vector)

    obslon = Float64[]
    obslat = Float64[]
    obsdepth = Float64[]
    obsdates = DateTime[]
    obsval1 = Float64[]
    obsval2 = Float64[]

    for datafile in datafilelist
        lon, lat, dates, depth, T, S, depth_level_thickness = read_TS_data(datafile)
        append!(obslon, lon)
        append!(obslat, lat)
        append!(obsdepth, depth)
        append!(obsdates, dates)
        append!(obsval1, T)
        append!(obsval2, S)
    end

    return obslon::Vector{Float64}, obslat::Vector{Float64}, obsdates::Vector{DateTime}, 
    obsdepth::Vector{Float64}, obsval1::Vector{Float64}, obsval2::Vector{Float64}#, 
    #depth_level_thickness::Vector{Float64}

end


"""
    read_data(datafilelist)

Read all the observations from the list of files `datafilelist`, 
as obtained with the function `get_filelist`.

# Example
```julia-repl
julia> datafilelist = get_filelist("./data/ME4OHL/", 1999:2014)
julia> obslon, obslat, obsdepth, obsdates, T, S = read_data(datafilelist)
```
"""
function read_data(datafilelist::Vector{String})

    obslonall = Float64[]
    obslatall = Float64[]
    obsdepthall = Float64[]
    obsdatesall = DateTime[]
    Tall = Float64[]
    Sall = Float64[]

    for datafile in datafilelist
        lon, lat, dates, depth, T, S, dohc = read_profile(datafile)
        lon[lon .< 20.] .+= 360;
        obslon, obslat, obsdates, obsdepth, T, S = ME4OH.vectorize_obs(lon, lat, dates, depth, T, S);

        append!(obslonall, obslon)
        append!(obslatall, obslat)
        append!(obsdepthall, obsdepth)
        append!(obsdatesall, obsdates)
        append!(Tall, T)
        append!(Sall, S)
    end
    
    return obslonall::Vector{Float64}, obslatall::Vector{Float64}, obsdepthall::Vector{Float64},
        obsdatesall::Vector{DateTime}, Tall::Vector{Float64}, Sall::Vector{Float64}
end

"""
    read_data(datafilelist)

Read all the observations from the list of files `datafilelist`, 
as obtained with the function `get_filelist`.

# Example
```julia-repl
julia> datafilelist = get_filelist("./data/ME4OHL/", 1999:2014)
julia> obslon, obslat, obsdepth, obsdates, dohc = read_data_dohc(datafilelist)
```
"""
function read_data_dohc(datafilelist::Vector{String})
    obslonall = Float32[]
    obslatall = Float32[]
    obsdepthall = Float32[]
    obsdatesall = DateTime[]
    obsvalall = Float32[]

    for datafile in datafilelist
        lon, lat, dates, vertical_levels, T, S, dohc, dohc_mask, depthbounds = read_profile(datafile)
        lon[lon .< 20.] .+= 360;
        obslon, obslat, obsdepth, obsdates, obsval = vectorize_dohc(lon, lat, 
        depthbounds, dates, dohc, dohc_mask)

        append!(obslonall, obslon)
        append!(obslatall, obslat)
        append!(obsdepthall, obsdepth)
        append!(obsdatesall, obsdates)
        append!(obsvalall, obsval)
    end

    return obslonall::Vector{Float32}, obslatall::Vector{Float32}, obsdepthall::Vector{Float32},
        obsdatesall::Vector{DateTime}, obsvalall::Vector{Float32}
end

"""
    vectorize_obs(lon, lat, dates, depth, T, S)

Transform the observations (list of profiles) into vectors.
This implies repeating the coordinates vectors.

# Example
```julia-repl 
julia> obslon, obslat, obsdates, obsdepth, T, S = vectorize_obs(lon, lat, dates, depth, T, S)
```
"""
function vectorize_obs(lon, lat, dates, depth, T, S)

    nprofiles = length(lon)
    nlevels = length(depth)
    
    obslon = repeat(lon, inner=nlevels)
    obslat = repeat(lat, inner=nlevels)
    obsdepth = repeat(depth, outer=nprofiles)
    obsdates = repeat(dates, inner=nlevels);

    goodT = findall(.!isnan.(T[:]));
    obslon = obslon[goodT]
    obslat = obslat[goodT]
    obsdepth = obsdepth[goodT]
    obsdates = obsdates[goodT]
    T = T[goodT]
    S = S[goodT]
    return obslon::Vector{Float64}, obslat::Vector{Float64}, obsdates::Vector{DateTime}, 
    obsdepth::Vector{Float64}, T::Vector{Float64}, S::Vector{Float64}
end

"""
    vectorize_dohc(lon, lat, bounds, dates, dohc, dohc_mask)

Generate vectors of coordinates and observations for the `dohc` variable (defined on 3 levels).
The inputs are obtained with the function `read_profile`.

# Example
```julia-repl
julia> lon, lat, dates, vertical_levels, T, S, dohc, dohc_mask, bounds = read_profile(datafile)
julia> obslon, obslat, obsdepth, obsdates, obsval = vectorize_dohc(lon, lat, bounds, dates, dohc, dohc_mask);
```
"""
function vectorize_dohc(lon::Vector{Float32}, lat::Vector{Float32}, bounds::Matrix{Float32}, 
    dates::Vector{DateTime}, dohc::Matrix{Float32}, dohc_mask::AbstractArray{Bool})

nprofiles = length(lon)
ndepth = length(bounds[:,2])
ntotal = nprofiles * ndepth
lonrep = transpose(repeat(lon, outer=[1,ndepth]));
latrep = transpose(repeat(lat, outer=[1,ndepth]));
datesrep = Array{DateTime}(undef, ndepth, nprofiles)

for iii = 1:3
    datesrep[iii,:] = dates
end
depthrep = repeat(bounds[:,2], outer=[1, nprofiles])

# Find good values 
goodvalues = .!isnan.(dohc) .& dohc_mask
@info("Found $(sum(goodvalues)) out of $(ntotal) values")

obslon = lonrep[goodvalues]
obslat = latrep[goodvalues]
obsdates = datesrep[goodvalues]
obsdepth = depthrep[goodvalues]
obsval = dohc[goodvalues]

return obslon::Vector{Float32}, obslat::Vector{Float32}, obsdepth::Vector{Float32},
    obsdates::Vector{DateTime}, obsval::Vector{Float32}
end

"""
    me4oh_dohc(T, z)

Compute the ocean heat content density using the temperature field `T` and
the depth levels `z`.

# Example
```julia-repl 
julia> dohc = me4oh_dohc(T, Z)
```
"""
function me4oh_dohc(T, z)

    rho0 = 1030    # kg/m^3
    Cp0 = 3989.244 # J/(kg K)
    T0 = 273.15    # Celsius -> Kelvin
    tera = 10^12   # 
    scale = Cp0 * rho0/tera

    z = diff(z) * scale
    dohc = cumsum((T + T0) .* dz) # Km(kg/m^3)(J/kg K)/10^12 = TJ/m^2
    return dohc
end

"""
    get_year(datafile)

Return the year corresponding to the data file.     
Useful to list files for a given time interval.

# Example
```julia-repl
julia> get_year("ofam3-jra55.all.EN.4.1.1.f.profiles.g10.197901.update.nc)
```
"""
function get_year(datafile::AbstractString)
    NCDataset(datafile, "r") do ds
        theyear = Int64(ds["en4_ymd"][1,1])
        return theyear
    end
end

"""
    make_fname(timeperiod, depthlayer, experiment; product="DIVAnd")

Create the output netCDF file name based on the time period, the depth layer and the experiment

# Example
```julia-repl 
julia> make_fname(2005:2014, [286.6, 685.9], "A")
```
"""
function make_fname(timeperiod::UnitRange{Int64}, depthlayer::Vector{Float64}, experiment; product::String="DIVAnd")
    fname = "OHC_$(timeperiod[1])_$(timeperiod[end])_lev$(depthlayer[1])_$(depthlayer[end])_exp$(experiment)_$(product).nc"
    return fname
end

"""
    get_timegrid(timeperiod)

Create a time vector in the interval define by `timeperiod`, with a monthly resolution.

# Example
```julia-repl 
julia> datesr = get_timegrid(timeperiod
```
"""
function get_timegrid(timeperiod::UnitRange{Int64})
    dayref = 1 # or shoule be 15?
    timegrid = collect(DateTime(timeperiod[1], 1, dayref):Dates.Month(1):DateTime(timeperiod[end], 12, dayref))
    return timegrid::Vector{DateTime}
end

"""
    datetime2days(timeperiod, dateref = dateref)

Convert the DataTime vector to a vector of Float64,
to be used as an input in the netCDF creation.

# Example
```julia-repl 
julia> daygrid = datetime2days(timeperiod)
```
"""
function datetime2days(timeperiod::UnitRange{Int64}; dateref = DateTime(1900, 1, 1))
    factor = 1000 * 3600 * 24
    timegrid = get_timegrid(timeperiod);
    daygrid = [(tt .- dateref).value / (factor) for tt in timegrid]
    return daygrid::Vector{Float64}
end

"""
    get_profile_number(datafilelist)

Compute the number of profiles in each file contained in `datafilelist`

# Example
```julia-repl 
julia> nprofiles = get_profile_number(datafilelist)
```
"""
function get_profile_number(datafilelist::Vector{String})
    nprofiles = zeros(Int64, length(datafilelist))
    for (ntime, datafile) in enumerate(datafilelist)
        
        # Read observations
        lon, lat, dates, depth, T, S, dohc = read_profile(datafile)
        nprofiles[ntime] = length(lon)
    end
    return nprofiles::Vector{Int64}
end

"""
    get_clim_bounds(timeperiod)

Compute the climatological bounds from the selected time period

# Example
```julia-repl 
julia> get_clim_bounds(1974:2015)
2×12 Matrix{DateTime}:
 1974-01-01T00:00:00  1974-02-01T00:00:00  1974-03-01T00:00:00  …  1974-10-01T00:00:00  1974-11-01T00:00:00  1974-12-01T00:00:00
 2015-01-31T00:00:00  2015-02-28T00:00:00  2015-03-31T00:00:00     2015-10-31T00:00:00  2015-11-30T00:00:00  2015-12-31T00:00:00
```
"""
function get_clim_bounds(timeperiod::UnitRange{Int64})
    climatology_bounds = Matrix{Dates.DateTime}(undef, 2, 12)

    for mm = 1:12
        climatology_bounds[1,mm] = DateTime(timeperiod[1], mm, 1)
        climatology_bounds[2,mm] = DateTime(timeperiod[end], mm, Dates.daysinmonth(timeperiod[end], mm))
    end

    return climatology_bounds::Matrix{Dates.DateTime}
end

"""
    create_netcdf_results(fname, longrid, latgrid, timegrid, valex, title)

Create the netCDF file with the spatial (defined by `longrid` and `latgrid`) and temporal grid (defined by `timegrid`).

According to the documentation `ME4OH_protocol_FINAL`: 
> One NetCDF file for each depth layer

# Example
```julia-repl 
julia> create_netcdf_results("results.nc", "dohc", -180:0.5:180., -75.:0.5:75., 2005:2014)
```
"""
function create_netcdf_results(fname::AbstractString, varname::String, longrid::StepRangeLen, latgrid::StepRangeLen, timeperiod::UnitRange{Int64}; 
    valex=-999, title::AbstractString="DIVAnd interpolated field")
    
    daygrid = datetime2days(timeperiod);
    globalattribs = OrderedDict(
        "creation_date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "title" => title,
        "institute" => "GHER, FOCUS, University of Liège",
        "Tool" => "DIVAnd",
        "Tool version" => "2.7.11",
        "Julia version" => "1.11.0-rc2",
        "Author" => "Charles Troupin",
        "Author email" => "ctroupin@uliege.be",
        "Author orcID" => "0000-0002-0265-1021",  
    )

    NCDataset(fname, "c", attrib = globalattribs) do ds

    # Dimensions
    ds.dim["lat"] = length(latgrid)
    ds.dim["lon"] = length(longrid)
    ds.dim["time"] = Inf

    # Declare variables
    defVar(ds, "lat", latgrid, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "actual_range"              => [minimum(latgrid), maximum(latgrid)],
        "long_name"                 => "Latitude",
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))

    defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "actual_range"              => [minimum(longrid), maximum(longrid)],
        "long_name"                 => "Longitude",
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    defVar(ds, "time", daygrid, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "actual_range"              => [Dates.format(Date(timeperiod[1],1,1), dateformat"yyyy-mm-dd"), Dates.format(Date(timeperiod[2],12,31), dateformat"yyyy-mm-dd")],
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    # Create variable storing the gridded field
    # Units and names change according to the variable 
    if varname == "dohc"
        defVar(ds, varname, Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "TJ/m^2",
            "short_name"                => "ocean_heat_content_density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

        defVar(ds, varname * "_error", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "1",
            "short_name"                => "ocean_heat_content_density_error",
            "long_name"                 => "Relative error on the heat content density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

    elseif varname == "adohc"
        defVar(ds, varname, Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "GJ/m^2",
            "short_name"                => "ocean_heat_content_density_anomaly",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

        defVar(ds, varname * "_error", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "1",
            "short_name"                => "ocean_heat_content_density_anomaly_error",
            "long_name"                 => "Relative error on the heat content density anomaly",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

    elseif varname == "dadohc"
        defVar(ds, varname, Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "GJ/m^2",
            "short_name"                => "detrended_ocean_heat_content_density_anomaly",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

        defVar(ds, varname * "_error", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "GJ/m^2",
            "short_name"                => "detrended_ocean_heat_content_density_anomaly_error",
            "long_name"                 => "Relative error on the heat content density detrended anomaly",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))

    else 
        @error("Variable $(varname) does not exist")
    end
    

    return nothing
    end;
end


"""
    create_netcdf_climatology(fname, longrid, latgrid, timegrid)

Create the netCDF file with the spatial (defined by `longrid` and `latgrid`) that will store the monthly climatology.

# Example
```julia-repl 
julia> create_netcdf_climatology("results.nc", "dohc", -180:0.5:180., -75.:0.5:75., 2005:2014)
```
"""
function create_netcdf_climatology(fname::AbstractString, varname::String, 
    longrid::StepRangeLen, latgrid::StepRangeLen, timeperiod::UnitRange{Int64}; 
    valex=-999, title::AbstractString="DIVAnd interpolated field")
    
    meanyear = Int64(floor(0.5 * (timeperiod[1] + timeperiod[end])))
    daygrid = datetime2days(meanyear:meanyear);

    globalattribs = OrderedDict(
        "creation_date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "title" => title,
        "institute" => "GHER, FOCUS, University of Liège",
        "Tool" => "DIVAnd",
        "Tool version" => "2.7.11",
        "Julia version" => "1.11.0-rc2",
        "Author" => "Charles Troupin",
        "Author email" => "ctroupin@uliege.be",
        "Author orcID" => "0000-0002-0265-1021",  
    )

    NCDataset(fname, "c", attrib = globalattribs) do ds

    # Dimensions
    ds.dim["level_bound"] = 3
    ds.dim["nv"] = 2
    ds.dim["lat"] = length(latgrid)
    ds.dim["lon"] = length(longrid)
    ds.dim["time"] = 12

    # Declare variables
    defVar(ds, "lat", latgrid, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "actual_range"              => [minimum(latgrid), maximum(latgrid)],
        "long_name"                 => "Latitude",
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))

    defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "actual_range"              => [minimum(longrid), maximum(longrid)],
        "long_name"                 => "Longitude",
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    defVar(ds, "time", daygrid, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    climatology_bounds = get_clim_bounds(timeperiod)

    defVar(ds,"climatology_bounds", climatology_bounds, ("nv", "time"), attrib = OrderedDict(
        "units"                     => "days since 1900-01-01 00:00:00",
    ))


    # Create variable storing the gridded field
    defVar(ds, "dohc", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "TJ/m^2",
            "short_name"                => "ocean_heat_content_density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
    ))

    defVar(ds, "dohc_error", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "1",
            "short_name"                => "ocean_heat_content_density",
            "long_name"                 => "Relative error on the heat content density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
    ))
    

    return nothing
    end;
end

"""
    create_netcdf_climatology4D(fname, longrid, latgrid, depthgrid, timegrid)

Create the netCDF file with the spatial (defined by `longrid`, `latgrid`, `depthgrid`) that will store 
the monthly climatology.

# Example
```julia-repl 
julia> create_netcdf_climatology("results.nc", -180:0.5:180., -75.:0.5:75., [2.5, 5., 10., 25.], 2005:2014)
```
"""
function create_netcdf_climatology4D(fname::AbstractString, longrid::StepRangeLen, latgrid::StepRangeLen, 
    depthgrid::Vector{Float64}, timeperiod::UnitRange{Int64}; 
    valex=-999, title::AbstractString="DIVAnd temperature interpolated field")
    
    meanyear = Int64(floor(0.5 * (timeperiod[1] + timeperiod[end])))
    daygrid = datetime2days(meanyear:meanyear);

    globalattribs = OrderedDict(
        "creation_date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "title" => title,
        "institute" => "GHER, FOCUS, University of Liège",
        "Tool" => "DIVAnd",
        "Tool version" => "2.7.11",
        "Julia version" => "1.11.0-rc2",
        "Author" => "Charles Troupin",
        "Author email" => "ctroupin@uliege.be",
        "Author orcID" => "0000-0002-0265-1021",  
    )

    NCDataset(fname, "c", attrib = globalattribs) do ds

    # Dimensions
    ds.dim["depth"] = length(depthgrid)
    ds.dim["nv"] = 2
    ds.dim["lat"] = length(latgrid)
    ds.dim["lon"] = length(longrid)
    ds.dim["time"] = 12

    # Declare variables
    defVar(ds, "lat", latgrid, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "actual_range"              => [minimum(latgrid), maximum(latgrid)],
        "long_name"                 => "Latitude",
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))

    defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "actual_range"              => [minimum(longrid), maximum(longrid)],
        "long_name"                 => "Longitude",
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    defVar(ds, "time", daygrid, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))


    defVar(ds, "depth", depthgrid, ("depth",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Depth",
        "axis"                      => "Z",
        "long_name"                 => "Depth of measurement",
        "standard_name"             => "depth",
        "units"                     => "m",
    ))

    climatology_bounds = get_clim_bounds(timeperiod)

    defVar(ds,"climatology_bounds", climatology_bounds, ("nv", "time"), attrib = OrderedDict(
        "units"                     => "days since 1900-01-01 00:00:00",
    ))

    # Create variable storing the gridded field
    defVar(ds, "temperature", Float64, ("lon", "lat", "depth", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "degree_Celsius",
            "short_name"                => "sea_water_temperature",
            "standard_name"             => "sea_water_temperature"
    ))

    defVar(ds, "temperature_error", Float64, ("lon", "lat", "depth", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "1",
            "short_name"                => "sea water temperature error",
            "long_name"                 => "Relative error on the sea water temperature",
            "standard_name"             => "sea_water_temperature"
    ))
    

    return nothing
    end;
end


"""
    create_netcdf_product4D(fname, longrid, latgrid, depthgrid, timegrid)

Create the netCDF file with the spatial (defined by `longrid`, `latgrid`, `depthgrid`) that will store 
the monthly climatology.

# Example
```julia-repl 
julia> create_netcdf_product4D("results.nc", -180:0.5:180., -75.:0.5:75., [2.5, 5., 10., 25.], 2005:2014)
```
"""
function create_netcdf_product4D(fname::AbstractString, longrid::StepRangeLen, latgrid::StepRangeLen, 
    depthgrid::Vector{Float64}, timeperiod::UnitRange{Int64}; 
    valex=-999, title::AbstractString="DIVAnd temperature interpolated field")
    
    daygrid = datetime2days(timeperiod);

    globalattribs = OrderedDict(
        "creation_date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "title" => title,
        "institute" => "GHER, FOCUS, University of Liège",
        "Tool" => "DIVAnd",
        "Tool version" => "2.7.11",
        "Julia version" => "1.11.0-rc2",
        "Author" => "Charles Troupin",
        "Author email" => "ctroupin@uliege.be",
        "Author orcID" => "0000-0002-0265-1021",  
    )

    NCDataset(fname, "c", attrib = globalattribs) do ds

    # Dimensions
    ds.dim["depth"] = length(depthgrid)
    ds.dim["nv"] = 2
    ds.dim["lat"] = length(latgrid)
    ds.dim["lon"] = length(longrid)
    ds.dim["time"] = Inf

    # Declare variables
    defVar(ds, "lat", latgrid, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "actual_range"              => [minimum(latgrid), maximum(latgrid)],
        "long_name"                 => "Latitude",
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))

    defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "actual_range"              => [minimum(longrid), maximum(longrid)],
        "long_name"                 => "Longitude",
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    defVar(ds, "time", daygrid, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))


    defVar(ds, "depth", depthgrid, ("depth",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Depth",
        "axis"                      => "Z",
        "long_name"                 => "Depth of measurement",
        "standard_name"             => "depth",
        "units"                     => "m",
    ))

    climatology_bounds = get_clim_bounds(timeperiod)

    defVar(ds,"climatology_bounds", climatology_bounds, ("nv", "time"), attrib = OrderedDict(
        "units"                     => "days since 1900-01-01 00:00:00",
    ))

    # Create variable storing the gridded field
    defVar(ds, "temperature", Float64, ("lon", "lat", "depth", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "degree_Celsius",
            "short_name"                => "sea_water_temperature",
            "standard_name"             => "sea_water_temperature"
    ))

    defVar(ds, "temperature_error", Float64, ("lon", "lat", "depth", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "1",
            "short_name"                => "sea water temperature error",
            "long_name"                 => "Relative error on the sea water temperature",
            "standard_name"             => "sea_water_temperature"
    ))
    
    return nothing
    end;
end


"""
    add_residuals(outputfile::AbstractString, datafilelist)

Add the residual variables (not the values into an existing netCDF file, 
as created by the function `create_netcdf_results`.

# Example
```julia-repl 
julia> 
```
"""
function add_residuals(outputfile::AbstractString, datafilelist::Vector{String}; valex=-999.)

    # Need to know the number of profiles for each time period (i.e. each input file)
    nprofiles = get_profile_number(datafilelist)
    nprofmax = maximum(nprofiles)

    NCDataset(outputfile, "a") do ds
        
        ds.dim["nprofiles"] = nprofmax

        defVar(ds, "obslat", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "axis"                      => "Y",
            "long_name"                 => "Latitude of the observations",
            "standard_name"             => "latitude",
            "units"                     => "degrees_north",
        ))
        
        defVar(ds, "obslon", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "axis"                      => "Y",
            "long_name"                 => "Longitude of the observations",
            "standard_name"             => "latitude",
            "units"                     => "degrees_north",
        ))

        defVar(ds, "obstime", Float64, ("nprofiles", "time",), attrib = OrderedDict(
            "_CoordinateAxisType"       => "Time",
            "axis"                      => "T",
            "calendar"                  => "Gregorian",
            "long_name"                 => "Time of the observations",
            "standard_name"             => "time",
            "time_origin"               => "01-JAN-1900 00:00:00",
            "units"                     => "days since 1900-01-01T00:00:00Z",
        ))
        
        defVar(ds, "dohc_residuals", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "TJ/m^2",
            "short_name"                => "ocean_heat_content_density_residuals",
            "long_name"                 => "Residuals of the ocean heat content density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))
    end

    return nothing

end

end