# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("u")

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

  us=ncvar_get(nc,"uxy_s")
  us.egrip=us[egrip.i,egrip.j,]
  us.ngrip=us[ngrip.i,ngrip.j,]
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)
  nom.egrip[i,]=ss.norm[i]*us.egrip
  nom.ngrip[i,]=ss.norm[i]*us.ngrip
  print(i)
}

mean.us.egrip=c()
for(t in 1:length(time)){
  mean.us.egrip[t]=sum(nom.egrip[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
mean.us.ngrip=c()
for(t in 1:length(time)){
  mean.us.ngrip[t]=sum(nom.ngrip[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-us-egrip.txt",sep="")
write.table(mean.us.egrip, file.out, quote = F, sep = " ", row.names=F)
file.out=paste(work.fldr,"/scripts/score-weighting/mean-us-ngrip.txt",sep="")
write.table(mean.us.ngrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom.egrip=nom.ngrip=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  us=ncvar_get(nc,"uxy_s")
  us.egrip=us[egrip.i,egrip.j,]
  us.ngrip=us[ngrip.i,ngrip.j,]
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
  nom.egrip[i,]=ss.norm[i]*(us.egrip-mean.us.egrip)^2
  nom.ngrip[i,]=ss.norm[i]*(us.ngrip-mean.us.ngrip)^2
  print(i)
}

sd.us.egrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.us.egrip[t]=sqrt(sum(nom.egrip[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-us-egrip.txt",sep="")
write.table(sd.us.egrip, file.out, quote = F, sep = " ", row.names=F)

sd.us.ngrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.us.ngrip[t]=sqrt(sum(nom.ngrip[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-us-ngrip.txt",sep="")
write.table(sd.us.ngrip, file.out, quote = F, sep = " ", row.names=F)

#######################################################################################################
#  U NEGIS
# score-weighted mean
nom.negis=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  us=ncvar_get(nc,"uxy_s")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  us.negis=c()
  for(t in 1:length(time)){
    us.negis[t]=mean(us[,,t][basins>8.5 & basins<9.5])
  }
  pd=length(time)
  nc_close(nc)
  nom.negis[i,]=ss.norm[i]*us.negis
  print(i)
}

mean.us.negis=c()
for(t in 1:length(time)){
  mean.us.negis[t]=sum(nom.negis[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-us-negis.txt",sep="")
write.table(mean.us.negis, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom.negis=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="") 
  nc = nc_open(file)
  us=ncvar_get(nc,"uxy_s")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  us.negis=c()
  for(t in 1:length(time)){
    us.negis[t]=mean(us[,,t][basins>8.5 & basins<9.5])
  }
  pd=length(time)
  nc_close(nc)
  nom.negis[i,]=ss.norm[i]*(us.negis-mean.us.negis)^2
  print(i)
}

sd.us.negis=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.us.negis[t]=sqrt(sum(nom.negis[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-us-negis.txt",sep="")
write.table(sd.us.negis, file.out, quote = F, sep = " ", row.names=F)


############ u averaged in NGRIP

# score-weighted mean
nom.ngrip=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  us=ncvar_get(nc,"uxy_bar")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  us.ngrip=us[ngrip.i,ngrip.j,]
  pd=length(time)
  nc_close(nc)
  nom.ngrip[i,]=ss.norm[i]*us.ngrip
  print(i)
}

mean.us.ngrip=c()
for(t in 1:length(time)){
  mean.us.ngrip[t]=sum(nom.ngrip[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-um-ngrip.txt",sep="")
write.table(mean.us.ngrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom.ngrip=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  us=ncvar_get(nc,"uxy_bar")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  us.ngrip=us[ngrip.i,ngrip.j,]
  pd=length(time)
  nc_close(nc)
  nom.ngrip[i,]=ss.norm[i]*(us.ngrip-mean.us.ngrip)^2
  print(i)
}

sd.us.ngrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.us.ngrip[t]=sqrt(sum(nom.ngrip[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-um-ngrip.txt",sep="")
write.table(sd.us.ngrip, file.out, quote = F, sep = " ", row.names=F)






