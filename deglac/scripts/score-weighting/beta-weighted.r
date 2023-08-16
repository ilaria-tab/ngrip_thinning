# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("beta")

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

  beta=ncvar_get(nc,"beta")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  beta.egrip=beta[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*beta.egrip
  print(i)
}

mean.beta.egrip=c()
for(t in 1:length(time)){
  mean.beta.egrip[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-beta-egrip.txt",sep="")
write.table(mean.beta.egrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  beta=ncvar_get(nc,"beta")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  beta.egrip=beta[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(beta.egrip-mean.beta.egrip)^2
  print(i)
}

sd.beta.egrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.beta.egrip[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-beta-egrip.txt",sep="")
write.table(sd.beta.egrip, file.out, quote = F, sep = " ", row.names=F)


#########################################################################################################################
# score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  beta=ncvar_get(nc,"beta")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  beta.negis=c()
  for(t in 1:length(time)){
   beta.negis[t]=mean(beta[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*beta.negis
  print(i)
}

mean.beta.negis=c()
for(t in 1:length(time)){
  mean.beta.negis[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-beta-negis.txt",sep="")
write.table(mean.beta.negis, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  beta=ncvar_get(nc,"beta")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  for(t in 1:length(time)){
    beta.negis[t]=mean(beta[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(beta.negis-mean.beta.negis)^2
  print(i)
}

sd.beta.negis=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.beta.negis[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-beta-negis.txt",sep="")
write.table(sd.beta.negis, file.out, quote = F, sep = " ", row.names=F)



