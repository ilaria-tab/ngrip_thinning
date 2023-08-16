# rmse for pd

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: ds-ne")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par

# Step 1: Calulate rmse for median sim
#mp = median.par
#file = paste(out.fldr,"/itmc.",itmc[mp],".btq.",q[mp],".cfstrm.",cfstr[mp],".nffdlt.",nff[mp],".enhshr.",enh[mp],"/yelmo2D.nc", sep="")
#d.s = raster(file, varname="z_srf_pd_err", band=nbands)
#sq.s = d.s^2
#n = length(d.s[values(d.s)!="NA"])  # total domain
#rmse.med = sqrt(sum(as.vector(sq.s),na.rm=T)/n)
#sigma=rmse.med


# Step 2: Calulate misfits
rmse.h=rmse.s=c()
for(i in 1:nr){

  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",i-1,"/yelmo2D.nc",sep="")
 
  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.h[i]=NA
    rmse.s[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.h[i]=NA
      rmse.s[i]=NA
      next
    }
  }
  nc_close(nc)

  d.s = raster(file, varname="z_srf_pd_err", band=nbands)
  basins = raster(file, varname="basins")
  #d.s = d.s[basins==1.4 | basins==2.1 | basins==2.2  | basins==9.1 | basins==9.2]
  d.s = d.s[(basins>1.3 & basins<1.5) | (basins>2.0 & basins<2.3) | (basins>9.0 & basins<9.3)]
  sq.s = d.s^2
  n = length(d.s)  # total domain

  d.h = raster(file, varname="H_ice_pd_err", band=nbands)
  basins = raster(file, varname="basins")
  d.h = d.h[(basins>1.3 & basins<1.5) | (basins>2.0 & basins<2.3) | (basins>9.0 & basins<9.3)]
  sq.h = d.h^2
  n = length(d.h)

  # Misfit 
  #c = -1/(2*n)
  #err = sum(as.vector(sq.s)/sigma^2,na.rm=T)
  #mis = c * err
  #mis.s[i] = mis
  
  rmse.s[i] = sqrt(sum(sq.s)/n)
  rmse.h[i] = sqrt(sum(sq.h)/n)
  print(i)

}

#b=cbind(a,rmse.s)

## write results to a file
file.out=paste(work.fldr,"/rmse-h-ne.txt",sep="")
write.table(rmse.h, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-s-ne.txt",sep="")
write.table(rmse.s, file.out, quote = F, sep = " ", row.names=F)


