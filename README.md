# VICmodel

## Overview
VICmodel is an R implementation of the Variable Infiltration Capacity (VIC) macroscale distributed hydrologic model (Liang et al., 1994) originally developed by Xu Liang at the University of Washington ([UW](http://www.washington.edu/)), USA, and currently the model is maintained by the Computational Hydrology group ([UW Hydro](http://uw-hydro.github.io/)) in the [Department of Civil and Environmental Engineering](https://www.ce.washington.edu/) at UW. This R package is developed by Ruida Zhong at the Center for Water Resources and Environment, [Sun Yat-sen University](http://www.sysu.edu.cn/) (SYSU). This package is built to for the more convinient use for the R users and other users or researchers using windows platform.

The VIC model can simulate several land surface processes physically based on both water balance and energy balance, e.g. Evapotranspiration (on vegetation canopy, vegetation transpiration and soil evaporation), runoff (surface and underground), changes of soil moisture, soil ice and soil temperature of each soil layer, accumulation and melt of snow, sensible and latent heat flux between atmosphere and land surface, streamflow of the basin outlet (needed to be coupled with a runoff routing model), and many other variables. The landsurface parameters (about vegetation, soil, topography) and the timeseries of meteorological data (meteorological forcing, including precipitation, air temperature, incomming shortwave and longwave radiation, wind speed, air pressure and vapor pressure) are necessary inputs to run the VIC model.

For more information about VIC please see the [official documentation website of VIC](http://vic.readthedocs.io/en/master/).

## Dependencies

**R** >= 3.0.0
 
**R-packages:**

- Rcpp >= 0.12.0
- foreach

## Installation

You can install VICmodel from github with:

``` r
# Via devtools
require(devtools)
devtools::install_github("Sibada/VICmodel")
```

## References
Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges (1994), A simple hydrologically based model of land surface water and energy fluxes for general circulation models, _J. Geophys. Res_., **99**(D7), 14415-14428, [doi:10.1029/94JD00483](http://dx.doi.org/10.1029/94JD00483).

## Example

This is an example to run the VIC model using the sample inputs:

``` r
#Sample data, with 16 gridcells and 10 day hourly meteorological forcing inputs
data(STEHE)

forcing <- STEHE$forcing
soil <- STEHE$soil
veg <- STEHE$veg

# Options and settings of the VIC model
vic_param('start_year', 1949)
vic_param('start_month', 1)
vic_param('start_day', 1)
vic_param('end_year', 1949)
vic_param('end_month', 1)
vic_param('end_day', 10)
vic_param('step_per_day', 24)
vic_param('snow_step_per_day', 24)
vic_param('runoff_step_per_day', 24)

# Definition of model outputs (output variables, timescale, etc.)
out_info <- list(
  wb = list(timescale = 'hour', aggpar = 6,
            outvars = c('OUT_RUNOFF', 'OUT_BASEFLOW', 'OUT_SOIL_MOIST'),
            aggtypes = c('sum', 'sum', 'end')),
  eb = list(timescale = 'day', aggpar = 1,
            outvars = c('OUT_SWE', 'OUT_SOIL_TEMP'),
            aggtypes = c('avg', 'min'))
)

# Run the VIC model
outputs <- vic(forcing, soil, veg, output_info = out_info)

# Show information of model run and outputs.
print(outputs)
```
