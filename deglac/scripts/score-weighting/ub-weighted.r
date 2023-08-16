# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("ub")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# score-weighted mean
nom.egrip=nom.ngrip=array(0, dim=c(nr.best,nbands))
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

#nc = nc_open(file)
  ub=ncvar_get(nc,"uxy_b")
  ub.egrip=ub[egrip.i,egrip.j,]
  ub.ngrip=ub[ngrip.i,ngrip.j,]
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)
  nom.egrip[i,]=ss.norm[i]*ub.egrip
  nom.ngrip[i,]=ss.norm[i]*ub.ngrip
  print(i)
}

mean.ub.egrip=c()
for(t in 1:length(time)){
  mean.ub.egrip[t]=sum(nom.egrip[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
mean.ub.ngrip=c()
for(t in 1:length(time)){
  mean.ub.ngrip[t]=sum(nom.ngrip[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-ub-egrip.txt",sep="")
write.table(mean.ub.egrip, file.out, quote = F, sep = " ", row.names=F)
file.out=paste(work.fldr,"/scripts/score-weighting/mean-ub-ngrip.txt",sep="")
write.table(mean.ub.ngrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom.egrip=nom.ngrip=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  ub=ncvar_get(nc,"uxy_b")
  ub.egrip=ub[egrip.i,egrip.j,]
  ub.ngrip=ub[ngrip.i,ngrip.j,]
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  pd=length(time)
  nc_close(nc)
  nom.egrip[i,]=ss.norm[i]*(ub.egrip-mean.ub.egrip)^2
  nom.ngrip[i,]=ss.norm[i]*(ub.ngrip-mean.ub.ngrip)^2
  print(i)
}

sd.ub.egrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.ub.egrip[t]=sqrt(sum(nom.egrip[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-ub-egrip.txt",sep="")
write.table(sd.ub.egrip, file.out, quote = F, sep = " ", row.names=F)

sd.ub.ngrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.ub.ngrip[t]=sqrt(sum(nom.ngrip[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-ub-ngrip.txt",sep="")
write.table(sd.ub.ngrip, file.out, quote = F, sep = " ", row.names=F)

#######################################################################################################
#  U NEGIS
# score-weighted mean
nom.negis=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  ub=ncvar_get(nc,"uxy_b")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  ub.negis=c()
  for(t in 1:length(time)){
    ub.negis[t]=mean(ub[,,t][basins>8.5 & basins<9.5])
  }
  pd=length(time)
  nc_close(nc)
  nom.negis[i,]=ss.norm[i]*ub.negis
  print(i)
}

mean.ub.negis=c()
for(t in 1:length(time)){
  mean.ub.negis[t]=sum(nom.negis[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-ub-negis.txt",sep="")
write.table(mean.ub.negis, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom.negis=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  ub=ncvar_get(nc,"uxy_b")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  ub.negis=c()
  for(t in 1:length(time)){
    ub.negis[t]=mean(ub[,,t][basins>8.5 & basins<9.5])
  }
  pd=length(time)
  nc_close(nc)
  nom.negis[i,]=ss.norm[i]*(ub.negis-mean.ub.negis)^2
  print(i)
}

sd.ub.negis=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.ub.negis[t]=sqrt(sum(nom.negis[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-ub-negis.txt",sep="")
write.table(sd.ub.negis, file.out, quote = F, sep = " ", row.names=F)
