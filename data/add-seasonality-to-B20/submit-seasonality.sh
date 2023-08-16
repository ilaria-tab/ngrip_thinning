#!/bin/bash 
#SBATCH -p short
#SBATCH -J R
#SBATCH -N 1         # nodes number
#SBATCH -n 48        # cores (CPU) number
#SBATCH -t 5-00:00 
#SBATCH --mem=100000
#SBATCH -o r.out
#SBATCH -e r.err
module purge
module load gnu8 R/4.0.3
module load impi
module load netcdf/4.6.3

# Path for some libraries
#export NCDF4=/mnt/lustre/home/itabone/R/x86_64-redhat-linux-gnu-library/3.6/ncdf4
#export PATH=.:$NCVIEW_HOME/libs:$PATH

# Run the job
Rscript season-improved.r 

