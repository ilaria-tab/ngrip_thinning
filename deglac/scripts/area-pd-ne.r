# rmse for gl pd

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: area-pd-ne")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par

# Margin from TOPO 
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
mask.topo=raster(file, varname="mask")
mask.topo[mask.topo==0 | mask.topo==1]=0
#mask.topo=boundaries(mask.topo, "inner", classes=F)
#mask.topo[mask.topo==0]=NA  # Now mask.topo has the margins =1, while is NA anywhere else
mask.topo[mask.topo==2]=1
area.topo=mask.topo

# Mask Ellesmere island, E and NE
topo.regions="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_REGIONS.nc"
region=raster(topo.regions, varname="mask")  # Ellesmere 1.11
area.topo[region < 1.27 | region > 1.3]=0

# BAsins ==2 
basin.file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa.nc"
basins=raster(basin.file, varname="basin_sub")

# Calulate misfits
rmse.area.ne=c()

for(i in 1:nr){

  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")

  file=paste(out.fldr,"/",i-1,"/yelmo2D.nc",sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.area.ne[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.area.ne[i]=NA
      next
    #nc_close(nc)
    }
  }
  nc_close(nc)

  # create Yelmo grounded mask 
  fice=raster(file, varname="f_ice", band=nbands)
  mask=raster(file, varname="mask_bed", band=nbands)
  #fice[mask==5 | mask==0]=0 # if floating or ocean, set it to 0
  fice[mask==1 | mask==0]=0 # if land or ocean, set it to 0
  area.y=fice
  area.y[basins<1.3 | basins>2.3]=NA

  # mask for Ellesmere island, E, NE
  area.y[region < 1.27 | region > 1.3]=NA
  
  # calculate difference
  dif=area.y - area.topo   # dif=0 if they coincide (ocean, land or grounded), dif=1 if ice is only modelled, dif=-1 if ice is observed but not modelled
    
  dif.v=as.vector(dif)*8*8 # dif*8km*8km = dif in area (km^2)
  rmse.area.ne[i]=sqrt(mean(dif.v^2, na.rm=T)) 
  print(i)
}

## write results to a file
file.out=paste(work.fldr,"/rmse-area-pd-ne.txt",sep="")
write.table(rmse.area.ne, file.out, quote = F, sep = " ", row.names=F)
