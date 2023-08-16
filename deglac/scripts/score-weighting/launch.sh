#!/bin/bash

conda activate r_env

Rscript bedrock-weighted.r 
Rscript beta-weighted.r
Rscript gl-weighted.r                
Rscript hice-weighted.r
Rscript neff-weighted.r
Rscript smb-weighted.r
Rscript taub-weighted.r
Rscript thinning-weighted.r
Rscript u-weighted.r
Rscript ub-weighted.r

