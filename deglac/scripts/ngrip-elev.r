# NGRIP elevation change

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: ngrip-elev")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}

# Step 2: Calulate misfits
rmse.ngrip=c()
for (i in 1:nr){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
   file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.ngrip[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.ngrip[i]=NA
      next
    }
  }
  #nc_close(nc)

  dz.srf=ncvar_get(nc, "z_srf_pd_err")
  time.y=ncvar_get(nc,"time")
  #time.y=time.y-1950
  dz.ngrip=dz.srf[ngrip.i, ngrip.j,]
  #lines(time.y, z.ngrip, col="blue", lty=1, lwd=1)

  # RMSE
  t.i=which(time.y==-11800)
  t.f=length(time.y)
  #t.f=which(time.y==-4000)
  mod.elev=dz.ngrip[t.i:t.f]
  time.y=time.y[t.i:t.f]

  obs=c() 
  for(tt in 1:length(time.y)){
    #ii=which(t.ds.v09==time.y[tt])
    ii=which(abs(time.y[tt]-t.ds.v09)==min(abs(time.y[tt]-t.ds.v09)))
    obs[tt]=ds.ngrip.v09[ii]
  }

  d = mod.elev - obs
  sq.k = d^2
  n = length(d)

  # RMSE
  rmse.ngrip[i]=sqrt(sum(as.vector(sq.k),na.rm=T)/n)
  print(i)
  nc_close(nc)
}

#dev.off()

## write results to a file
file.out=paste(work.fldr,"/rmse-ngrip.txt",sep="")
write.table(rmse.ngrip, file.out, quote = F, sep = " ", row.names=F)

