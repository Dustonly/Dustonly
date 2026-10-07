## Overview
This version of the offline dust emission model, Dustonly, now contains changes done on the treatment of the mineralogical dataset. 

### Model Description
The main characteristics of the repository remain the same, that being:
Dustonly is the standalone version of the dust emission scheme used in the aerosol transport model MUSCAT ([Heinold et al.,2007](http://doi.org/10.1029/2006jd007443); [Wolke et al.,2012](http://doi.org/10.1016/j.atmosenv.2012.02.085)). The model is based on the dust emission model by [Tegen et al. (2002)](http://doi.org/10.1029/2001jd000963). It uses the emission scheme by [Marticorena and Bergametti (1995)](http://doi.org/10.1029/95jd00690).

### Installation
The model can be compiled with gfortran using the `make` comment. 
May include NetCDF in `$LDFLAGS` and `$CPPFLAGS`.
Otherwise, the NetCDF path can be specified in the Makefile.

### Usage
To run Dustonly a namelist called `INPUT` and some data fields are necessary.
A minimal setup can be found [here](https://doi.org/10.5281/zenodo.8320600).


## What's Changed
v1.0 with mineralogy included, mimicked the soil size distribution of the minerals described in the GMINER dataset ([Nickovic et al.,2012](https://doi.org/10.5194/acp-12-845-2012)) to the aerosol size distribution
v2.0 changes that by redistributing minerals to other sizes based on measurements and the Brittle Fragmentation Theory ([Kok et al.,2012](https://doi.org/10.5194/acp-12-845-2012); [Perlwitz et al., 2015a](https://doi.org/10.5194/acp-15-11593-2015); [Gonçalves et al., 2023]( https://doi.org/10.5194/acp-23-8623-2023)) for increased accuracy. 

* For a physical description of the changes, challenges, and validation, a preprint is avaliable at: [Gómez Maqueo Anaya et al., 2026](https://doi.org/10.5194/egusphere-2026-23)
* For details on the changes on the init_mineralmap subroutine visit https://github.com/Dustonly/Dustonly/pull/1

### Sub-grid FCOVER classes (`veg_classes`)
With `veg_classes = K > 0` (requires `vegdayFile`, `dust_scheme = 1`, `veg_scheme = 1, 2` or `3`) the daily
vegetation file holds the FCOVER distribution of each cell as written by
[dust-input-pipline](https://github.com/Dustonly/dust-input-pipline): `FCOVER_CLASS(time,class,rlat,rlon)` (mean
FCOVER of class k) and the class weights (`vegweight_name`, default `FCOVER_WEIGHT`); NaN (empty class) counts as
weight 0. The emission of a cell is `sum_k w_k * E(FCOVER_CLASS_k)`, i.e. the flux is computed once per class.
`FCOVER_WEIGHT` are fractions of the emitting area (sum 1, combined with the source map, `psrcType = 3`);
`vegweight_name = 'FCOVER_CELL_WEIGHT'` takes fractions of the whole cell (sum = emitting fraction, the rest of
the cell does not emit), to be used without source map. With `veg_classes = 0` (default) nothing changes.
