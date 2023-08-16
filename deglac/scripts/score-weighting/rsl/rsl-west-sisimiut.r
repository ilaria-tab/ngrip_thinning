# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

print("rsl-west-sisimiut")

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# Central West, SIsimiut (Bennike 2011, Long 2009)
loc=findij(66.953211, -53.697405)
t.min=c(-10271,-10587,-10169,-8524,-5283,-5881,-4817,-3578,-2343,-1345,-2323,-2305,-1281,-1486)
t.max=c(-10115,-10130,-8762,-8322,-4868,-5606,-4449,-3399,-2153,-1183,-2120,-2002,-1081,-1300)
rsl=c(129.9,94.9,98.4,39.8,18.4,13,7.68,0.69,-2.19,-4.22,-4.025,-4.025,-1.835,-3.865)
rsl.err=c(0.5,0.5,0.5,0.5,0.5,0.075,0.16,0.41,0.31,0.28,0.41,0.41,0.31,0.41)

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
file.out=paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-west-sisimiut.txt",sep="")
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
file.out=paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-west-sisimiut.txt",sep="")
write.table(sd.rsl.y, file.out, quote = F, sep = " ", row.names=F)

