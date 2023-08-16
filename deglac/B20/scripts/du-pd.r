library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: du-pd")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par

# Step 1: Calulate rmse for median sim
#mp = median.par
#file = paste(out.fldr,"/itmc.",itmc[mp],".btq.",q[mp],".cfstrm.",cfstr[mp],".nffdlt.",nff[mp],".enhshr.",enh[mp],"/yelmo2D.nc", sep="")
#d.h = raster(file, varname="H_ice_pd_err", band=nbands)
#sq.h = d.h^2
#n = length(d.h[values(d.h)!="NA"])  # total domain
#rmse.med = sqrt(sum(as.vector(sq.h),na.rm=T)/n)
#sigma=rmse.med

# Mask Ellesmere island
topo.regions="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_REGIONS.nc"
region=raster(topo.regions, varname="mask")  # Ellesmere 1.11

# Step 2: Calulate misfits
rmse.u=c()
for(i in 1:nr){

  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",i-1,"/yelmo2D.nc",sep="")
  
  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.u[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.u[i]=NA
      next
    }
  }
  nc_close(nc)

  d.u = raster(file, varname="uxy_s_pd_err", band=nbands)
  d.u[region < 1.12 & region > 1.10]=0 
  basins = raster(file, varname="basins")
  #d.u=d.u[basins==9]
  sq.u = d.u^2
  n = length(d.u[values(d.u)!="NA"])  # total domain

  # Misfit 
  #c = -1/(2*n)
  #err = sum(as.vector(sq.h)/sigma^2,na.rm=T)
  #mis = c * err
  #mis.h[i] = mis

  rmse.u[i] = sqrt(sum(as.vector(sq.u),na.rm=T)/n)
  print(i)
}

file.out=paste(work.fldr,"/rmse-u.txt",sep="")
write.table(rmse.u, file.out, quote = F, sep = " ", row.names=F)

