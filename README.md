
This repository contains the code for the creation of gridded fields of Density of Ocean Heat Content (`dohc`) in the frame of the ME4OH initiative, using the [`DIVAnd`](https://github.com/gher-uliege/DIVAnd.jl/) software tool. 

## DIVAnd 

`DIVAnd` stands for _Data-Interpolating Variational Analysis in n dimensions_ (Barth et al., 2014). The tool performs an n-dimensional variational analysis/gridding of arbitrarily located observations. It is a generalisation of [`DIVA`](https://github.com/gher-uliege/DIVA) (Troupin et al., 2012; Beckers et al., 2014), which was focused on two-dimensional interpolations only.

The main advantages of the method are: 
1. The computation time doesn't depend on the number of observations, make it possible to work on the World Ocean at several depth or for different time periods.
2. The interpolation takes into account physical boundaries, such as the coastline and the bathymetry, without having to apply a mask _a posteriori_.

### Implementation

### Usage

The method has been applied in several European projects such as SeaDataCloud, EMODnet Biology, EMODnet Chemistry and BlueCloud. It was also applied to produced regional and local climatologies (Belgacem et al., 2021; ). 



## Objectives

Different input files are provided and the goal is to interpolate the observations using different techniques (inclusing `DIVAnd`).

## Experiment A

The observations are interpolated layer by layer and time period by time period, i.e. we don't introduce any influence from one period to the next or from one layer to another.


## References

Barth, A., Beckers, J.-M., Troupin, C., Alvera-Azcárate, A. and Vandenbulcke, L.: divand-1.0: n-dimensional variational data analysis for ocean observations, Geoscientific Model Development, 7, 225–241, 2014.
DOI: [10.5194/gmd-7-225-2014](http://dx.doi.org/10.5194/gmd-7-225-2014) 

Belgacem, M., Schroeder, K., Barth, A., Troupin, C., Pavoni, B., Raimbault, P., Garcia, N., Borghini, M. and Chiggiato, J.: Climatological distribution of dissolved inorganic nutrients in the western Mediterranean Sea (1981–2017), _Earth System Science Data_, __13(12)__, 5915–5949, 2021.
DOI: [10.5194/essd-13-5915-2021](http://dx.doi.org/10.5194/essd-13-5915-2021) 

Shahzadi, K., Pinardi, N., Barth, A., Troupin, C., Lyubartsev, V. and Simoncelli, S.: A New Global Ocean Climatology, _Frontiers in Environmental Science_, __9__, 2021.
DOI: [10.3389/fenvs.2021.711363](http://dx.doi.org/10.3389/fenvs.2021.711363) 

