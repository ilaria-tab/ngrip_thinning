#!/bin/bash

conda activate r_env

#Rscript 79-margin.r 
Rscript area-lgm.r
Rscript area-pd-ne.r 
Rscript area-pd.r                
Rscript dh-pd.r
Rscript ds-ne.r
Rscript du-pd-ne.r
Rscript du-pd.r
Rscript ice-cores.r
Rscript larsen-morains.r
Rscript marine-cores.r
Rscript ng-icefree.r
Rscript ngrip-elev.r
Rscript score.r

