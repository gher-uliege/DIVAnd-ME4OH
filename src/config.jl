datadir = "/home/ctroupin/data/ME4OH/data/en4.1.1/1979-2014/full/update/"
outputdir = "../output/"

isdir(outputdir) ? @debug("Directory already exists") : mkpath(outputdir)

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

# Domain grid
longrid = 20.5:1.0:379.5 
latgrid = -89.5:1.0:89.5
bathyfile = "../data/gebco_30sec_16.nc"