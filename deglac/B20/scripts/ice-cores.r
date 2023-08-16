# Ice cores elev

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: ice-cores")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
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
ngrip.max=mis.cores=rmse.cores=c()
for (i in 1:nr){
  
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  
  file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.cores[i]=NA
    mis.cores[i]=NA
    ngrip.max[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.cores[i]=NA
      mis.cores[i]=NA
      ngrip.max[i]=NA
      next
    }
  }
 # nc_close(nc)

  dz.srf=ncvar_get(nc, "z_srf_pd_err")
  pd=nbands
  hol=70 # band when t=-4000
  ##dz.campcen.0=dz.srf[campcen.i, campcen.j,pd]
  ##dz.dye3.0=dz.srf[dye3.i, dye3.j,pd]
  ##dz.grip.0=dz.srf[grip.i, grip.j,pd]
  dz.ngrip.0=dz.srf[ngrip.i, ngrip.j,pd]
  dz.ngrip.hol=dz.srf[ngrip.i, ngrip.j,hol]
  ##dz0=c(dz.ngrip.0,dz.grip.0,dz.campcen.0,dz.dye3.0)   # m
  dz0=dz.ngrip.0  
  dzhol=dz.ngrip.hol

  sq.k = dz0^2
  n = length(dz0)
  #n = length(dzhol)  

  # max ngrip
  z.srf=ncvar_get(nc, "z_srf")
  z.ngrip=z.srf[ngrip.i,ngrip.j,]
  dz.ngrip=z.ngrip - z.ngrip[pd] 
  # ngrip.max[i]=max(dz.ngrip[11:nbands])
  ngrip.max[i]=max(dz.ngrip[11:nbands])

  # Misfit
  #c = -1/(2*n)
  #err = sum(as.vector(sq.k)/sigma^2,na.rm=T)
  #mis = c * err
  #mis.cores[i] = mis

  # RMSE
  rmse.cores[i]=sqrt(sum(as.vector(sq.k),na.rm=T)/n)
  print(i)
  nc_close(nc)
}

## write results to a file
#file.out=paste(work.fldr,"/mis-cores.txt",sep="")
#write.table(mis.cores, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-cores.txt",sep="")
write.table(rmse.cores, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/ngrip-max.txt",sep="")
write.table(ngrip.max, file.out, quote = F, sep = " ", row.names=F)

















