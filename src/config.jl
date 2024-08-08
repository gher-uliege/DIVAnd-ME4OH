using Downloads
using PyCall
mpl = pyimport("matplotlib")
mpl.style.use("./me4oh.mplstyle")


function make_dox_url(fileid::AbstractString)
    doxurl = "https://dox.uliege.be/index.php/s/$(fileid)/download"
    return doxurl::String
end

function download_check(datafile::AbstractString, datafileURL::AbstractString)

    if !isfile(datafile)
        Downloads.download(datafileURL, datafile)
    else
        @info "File already downloaded"
    end
end

# Directories
datadir = "/home/ctroupin/data/ME4OH/data/en4.1.1/1979-2014/full/update/"
datadirextra = "/home/ctroupin/data/ME4OH/data/en4.1.1/1979-2014/full/update/extra"

bathydir = "../data/"
datatestdir = "../data/test/"
outputdir = "../output/"
figdir = "../figures/"

isdir(outputdir) ? @debug("Directory already exists") : mkpath(outputdir)
isdir(datatestdir) ? @debug("Directory already exists") : mkpath(datatestdir)
isdir(figdir) ? @debug("Directory already exists") : mkpath(figdir)

# Files for the unit testing
datafiletest1 = joinpath(datatestdir, "ofam3-jra55.all.EN.4.1.1.f.profiles.g10.197901.update.nc")
datafiletest1URL = make_dox_url("rtQcFZDGszhhtfV")
datafiletest2 = joinpath(datatestdir, "ofam3-jra55.all.EN.4.1.1.f.profiles.g10.197902.update.nc")
datafiletest2URL = make_dox_url("87O2XkvhagaAQbD")
datafiletest3 = joinpath(datatestdir, "ofam3-jra55.all.EN.4.1.1.f.profiles.g10.197901.update.extra.nc")
datafiletest3URL = make_dox_url("47L09AV3RQQuaLB")

# GEBCO bathymetry
gebco04file = joinpath(bathydir, "gebco_30sec_4.nc")
gebco04fileURL = make_dox_url("RSwm4HPHImdZoQP")
gebco08file = joinpath(bathydir, "gebco_30sec_8.nc")
gebco08fileURL = make_dox_url("U0pqyXhcQrXjEUX")
gebco16file = joinpath(bathydir, "gebco_30sec_16.nc")
gebco16fileURL = make_dox_url("U0pqyXhcQrXjEUX")

# The following parameter values are taken from the shared document 
# "ME4OH_protocol_FINAL"

cp0 = 3989.244      # J/kg/K is heat capacity
rho0= 1030          # kg/m3 is density

# Time period
timeperiod1 = 1979:2014 # (whole period) 
timeperiod2 = 2005:2014 # (Argo period) 
timeperiod3 = 1993:2014 # (satellite period)

# Depth layers
depthlayer1 = [0, 286.6]     # m layer
depthlayer2 = [286.6, 685.9]  # m layer (i.e. starting from the lower bound of layer 286.6m, i.e. not including the layer associated with nominal model depth 286.6m)
depthlayer3 = [685.9, 1985.3]
depthlayers = [depthlayer1, depthlayer2, depthlayer3]

# Domain grid
longrid = 20.5:1.0:379.5 
latgrid = -89.5:1.0:89.5

bathyfile = gebco16file

download_check(bathyfile, gebco16fileURL);