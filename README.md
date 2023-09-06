# Dustonly
A Dust Emission Model

## Model Description 
Dustonly is the standalone version of the dust emission scheme used in the aerosol transport model MUSCAT ([Heinold et al.,2007](http://doi.org/10.1029/2006jd007443); [Wolke et al.,2012](http://doi.org/10.1016/j.atmosenv.2012.02.085)).
The model is based on the dust emission model by [Tegen et al. (2002)](http://doi.org/10.1029/2001jd000963). It uses the emission scheme by [Marticorena and Bergametti (1995)](http://doi.org/10.1029/95jd00690).

## Installation
The model can be compiled with gfortran using the `make` comment. 
May include NetCDF in `$LDFLAGS` and `$CPPFLAGS`.
Otherwise, the NetCDF path can be specified in the Makefile.

## Usage
To run Dustonly a namelist called `INPUT` and some data fields are necessary.
A minimal setup can be found [here](https://doi.org/10.5281/zenodo.8320600).

