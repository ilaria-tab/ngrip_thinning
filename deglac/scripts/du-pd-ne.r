library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: du-pd-ne")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par

# Mask Ellesmere island
topo.regions="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_REGIONS.nc"
region=raster(topo.regions, varname="mask")  # Ellesmere 1.11

# BAsins ==2 
basin.file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa-negis.nc"
basins=raster(basin.file, varname="basin")

# Step 2: Calulate misfits
rmse.u.ne=c()
for(i in 1:nr){

  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")

  file=paste(out.fldr,"/",i-1,"/yelmo2D.nc",sep="")  

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.u.ne[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.u.ne[i]=NA
      next
    #nc_close(nc)
    }
  }
  nc_close(nc)

  d.u = raster(file, varname="uxy_s_pd_err", band=nbands)
  #d.u[region < 1.12 & region > 1.10]=0 
  #d.u[basins!=9]=NA
  d.u[basins<1.3 | (basins>233 & basins<9.0)]=NA 
  sq.u = d.u^2
  n = length(which(is.na(as.vector(d.u))==F))

  rmse.u.ne[i] = sqrt(sum(as.vector(sq.u),na.rm=T)/n)
  print(i)
}

file.out=paste(work.fldr,"/rmse-u-ne.txt",sep="")
write.table(rmse.u.ne, file.out, quote = F, sep = " ", row.names=F)

