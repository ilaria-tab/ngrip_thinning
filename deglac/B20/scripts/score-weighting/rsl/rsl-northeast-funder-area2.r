# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

print("rsl-northeast-funder-area2.r")

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# North Funder et al., 2011 (Area 2 - 83°N - Cape Ole)
loc=findij(83.374692, -28.079069)
t=c(-7301,-10124,-9938,-9898,-9718,-9716,-9633,-9225,-8786,-8375,-8313,-8068,-6804,-6399,-4281)
t.err=c(123,162,229,267,169,171,111,190,195,163,112,98,104,97,129)
rsl=c(16,12,24,30,13,26,27,25,8,10,4,13,4,2,2)

# score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")

  # does the file exist?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    nom[i,]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      nom[i,]=NA
      next
    }
  } 
 
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  
  nom[i,]=ss.norm[i]*rsl.y
  print(i)
}

mean.rsl.y=c()
for(t in 1:length(time)){
  mean.rsl.y[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area2.txt",sep="")
write.table(mean.rsl.y, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")

  # does the file exist?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    nom[i,]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      nom[i,]=NA
      next
    }
  }

  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  
  nom[i,]=ss.norm[i]*(rsl.y-mean.rsl.y)^2
  print(i)
}

sd.rsl.y=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.rsl.y[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area2.txt",sep="")
write.table(sd.rsl.y, file.out, quote = F, sep = " ", row.names=F)

