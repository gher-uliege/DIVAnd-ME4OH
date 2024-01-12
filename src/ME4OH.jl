using NCDatasets
using DataStructures

"""
    read_profile(datafile)

Read the coordinates and the measurements from a data file

# Example
```julia-repl 
julia> lon, lat, dates, vertical_levels, T, S, dohc
```
"""
function read_profile(datafile::AbstractString)
    NCDataset(datafile, "r") do ds
        lon = ds["ts_lon"][:]
        lat = ds["ts_lat"][:]
        time = ds["en4_ymd"][:,:]
        dates = Dates.DateTime.(time[1,:], time[2,:], time[3,:])
        vertical_levels = ds["ts_z"][:]
        depth_level_thickness = ds["ts_dz"][:]
        T = ds["temp"][:,:]
        SSH = ds["eta_t"][:]
        SST = ds["sst"][:]
        S = ds["salt"][:,:]
        dohc = ds["dohc"][:,:]
        return lon::Vector{Float32}, lat::Vector{Float32}, dates::typeof(dates), 
        vertical_levels::Vector{Float32}, T::Matrix{Float32}, S::Matrix{Float32}, dohc::Matrix{Float32}
    end
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
    T = T[goodT]
    S = S[goodT]
    return obslon::Vector{Float32}, obslat::Vector{Float32}, obsdates::Vector{DateTime}, 
    obsdepth::Vector{Float32}, T::Vector{Float32}, S::Vector{Float32}
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
    T0 = 273.15    # Celcius -> Kelvin
    tera = 10^12   # 
    scale = Cp0 * rho0/tera

    z = diff(z) * scale
    dohc = cumsum((T + T0) .* dz) # Km(kg/m^3)(J/kg K)/10^12 = TJ/m^2
    return dohc
end

"""
    make_fname(timeperiod, depthlayer, experiment; product="DIVAnd")

Create the output netCDF file name based on the time period, the depthj layer and the experiment

# Example
```julia-repl 
julia> nprofiles = get_profile_number(datafilelist)
```
"""
function make_fname(timeperiod, depthlayer, experiment; product="DIVAnd")
    fname = "OHC_$(timeperiod[1])_$(timeperiod[end])_lev0_$(depthlayer[end])_exp$(experiment)_$(product).nc"
    return fname
end

"""
    get_timegrid(timeperiod)

Create a time vector in the interval define by `timeperiod`, with a monthly resolution.

# Example
```julia-repl 
julia> nprofiles = get_profile_number(datafilelist)
```
"""
function get_timegrid(timeperiod::UnitRange{Int64})
    dayref = 1 # or shoule be 15?
    timegrid = collect(DateTime(timeperiod[1], 1, dayref):Dates.Month(1):DateTime(timeperiod1[end], 12, dayref))
    return timegrid::Vector{DateTime}
end

"""
    datetime2days(timegrid, dateref = dateref)

Convert the DataTime vector to a vector of Float64,
to be used as an input in the netCDF creation.

# Example
```julia-repl 
julia> nprofiles = get_profile_number(datafilelist)
```
"""
function datetime2days(timegrid::Vector{DateTime}; dateref = DateTime(1900, 1, 1))
    timegrid1 = get_timegrid(timeperiod1);
    daygrid = [(tt .- dateref).value / (1000 * 3600 * 24) for tt in timegrid]
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
    create_netcdf_results(fname, longrid, latgrid, timegrid)

Create the netCDF file with the spatial (defined by `longrid` and `latgrid`) and temporal grid (defined by `timegrid`).

# Example
```julia-repl 
julia>
```
"""
function create_netcdf_results(fname::AbstractString, longrid, latgrid, timegrid::Vector{DateTime}; valex=-999)
    
    daygrid = datetime2days(timegrid);

    NCDataset(fname, "c", attrib = OrderedDict("title" => "DIVAnd interpolated field")) do ds

    # Dimensions
    ds.dim["lat"] = length(latgrid)
    ds.dim["lon"] = length(longrid)
    ds.dim["time"] = Inf

    # Declare variables
    nclat = defVar(ds, "lat", latgrid, ("lat",), attrib = OrderedDict(
        "axis"                      => "Y",
        "actual_range"              => [minimum(latgrid), maximum(latgrid)],
        "long_name"                 => "Latitude",
        "standard_name"             => "latitude",
        "units"                     => "degrees_north",
    ))

    nclon = defVar(ds, "lon", longrid, ("lon",), attrib = OrderedDict(
        "axis"                      => "X",
        "actual_range"              => [minimum(longrid), maximum(longrid)],
        "long_name"                 => "Longitude",
        "standard_name"             => "longitude",
        "units"                     => "degrees_east",
    ))

    nctime = defVar(ds, "time", daygrid, ("time",), attrib = OrderedDict(
        "_CoordinateAxisType"       => "Time",
        "actual_range"              => [Dates.format(timegrid[1], dateformat"yyyy-mm-dd"), Dates.format(timegrid[end], dateformat"yyyy-mm-dd")],
        "axis"                      => "T",
        "calendar"                  => "Gregorian",
        "long_name"                 => "Time of measurement",
        "standard_name"             => "time",
        "time_origin"               => "01-JAN-1900 00:00:00",
        "units"                     => "days since 1900-01-01T00:00:00Z",
    ))

    ncfield = defVar(ds, "dohc", Float64, ("lon", "lat", "time"), attrib = OrderedDict(
        "_FillValue"                => Float64(valex),
        "units"                     => "TJ/m^2",
		"short_name"                => "ocean_heat_content_density",
        "standard_name"                 => "sea_water_potential_temperature_expressed_as_heat_content"
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

        ncobslat = defVar(ds, "obslat", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "axis"                      => "Y",
            "long_name"                 => "Latitude of the observations",
            "standard_name"             => "latitude",
            "units"                     => "degrees_north",
        ))
        
        ncobslon = defVar(ds, "obslon", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "axis"                      => "Y",
            "long_name"                 => "Longitude of the observations",
            "standard_name"             => "latitude",
            "units"                     => "degrees_north",
        ))

        ncobstime = defVar(ds, "obstime", Float64, ("nprofiles", "time",), attrib = OrderedDict(
            "_CoordinateAxisType"       => "Time",
            "axis"                      => "T",
            "calendar"                  => "Gregorian",
            "long_name"                 => "Time of the observations",
            "standard_name"             => "time",
            "time_origin"               => "01-JAN-1900 00:00:00",
            "units"                     => "days since 1900-01-01T00:00:00Z",
        ))
        
        ncresiduals = defVar(ds, "dohc_residuals", Float64, ("nprofiles", "time"), attrib = OrderedDict(
            "_FillValue"                => Float64(valex),
            "units"                     => "TJ/m^2",
            "short_name"                => "ocean_heat_content_density_residuals",
            "long_name"                 => "Residuals of the ocean heat content density",
            "standard_name"             => "sea_water_potential_temperature_expressed_as_heat_content"
        ))
    end

    return nothing

end
