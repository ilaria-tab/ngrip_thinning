# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("taub")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# taub EGRIP
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

#  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  taub.egrip=taub[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*taub.egrip
  print(i)
}

mean.taub.egrip=c()
for(t in 1:length(time)){
  mean.taub.egrip[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-taub-egrip.txt",sep="")
write.table(mean.taub.egrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  taub.egrip=taub[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(taub.egrip-mean.taub.egrip)^2
  print(i)
}

sd.taub.egrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.taub.egrip[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-taub-egrip.txt",sep="")
write.table(sd.taub.egrip, file.out, quote = F, sep = " ", row.names=F)


#########################################################################################################################
# score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  taub.negis=c()
  for(t in 1:length(time)){
   taub.negis[t]=mean(taub[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*taub.negis
  print(i)
}

mean.taub.negis=c()
for(t in 1:length(time)){
  mean.taub.negis[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-taub-negis.txt",sep="")
write.table(mean.taub.negis, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="") 
  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  for(t in 1:length(time)){
    taub.negis[t]=mean(taub[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(taub.negis-mean.taub.negis)^2
  print(i)
}

sd.taub.negis=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.taub.negis[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-taub-negis.txt",sep="")
write.table(sd.taub.negis, file.out, quote = F, sep = " ", row.names=F)

#########################################################################################################################
# taub NGRIP score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  taub.ngrip=c()
  taub.ngrip=taub[ngrip.i,ngrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*taub.ngrip
  print(i)
}

mean.taub.ngrip=c()
for(t in 1:length(time)){
  mean.taub.ngrip[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-taub-ngrip.txt",sep="")
write.table(mean.taub.ngrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  taub=ncvar_get(nc,"taub")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  taub.ngrip=taub[ngrip.i,ngrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(taub.ngrip-mean.taub.ngrip)^2
  print(i)
}

sd.taub.ngrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.taub.ngrip[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-taub-ngrip.txt",sep="")
write.table(sd.taub.ngrip, file.out, quote = F, sep = " ", row.names=F)



