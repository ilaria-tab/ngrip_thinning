## README ##

Some info about the files contained in this folder.

tas_main....nc and pr_main_....nc : Original B20 data.

badgeley2020.f90 : source code to regrid B20 data to Yelmo grid. To be put in gridding/src/.

gridder_badgeley.f90 : compile this script to create the executable file to regrid B20 with gridding. To be put in gridding/.

badgeley2020.r : R script that adds ICE-5G surface elevations to B20 climate reconstructions. This script directly modifies the output file GRL-16KM_CLIM-RECON-B20.nc, therefore this should be run after gridding.

GRL-16KM_CLIM-RECON-B20.nc : B20 output file, already complete of surface elevations, ready to force Yelmo.

not_used/ : folder containing other data released along with B20. Not used here.
