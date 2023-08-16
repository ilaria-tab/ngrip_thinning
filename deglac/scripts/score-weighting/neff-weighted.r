# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("neff")

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

#nc = nc_open(file)
  n_eff=ncvar_get(nc,"N_eff")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  n_eff.egrip=n_eff[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*n_eff.egrip
  print(i)
}

mean.n_eff.egrip=c()
for(t in 1:length(time)){
  mean.n_eff.egrip[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-neff-egrip.txt",sep="")
write.table(mean.n_eff.egrip, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  n_eff=ncvar_get(nc,"N_eff")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  n_eff.egrip=n_eff[egrip.i,egrip.j,]
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(n_eff.egrip-mean.n_eff.egrip)^2
  print(i)
}

sd.n_eff.egrip=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.n_eff.egrip[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-n_eff-egrip.txt",sep="")
write.table(sd.n_eff.egrip, file.out, quote = F, sep = " ", row.names=F)


#########################################################################################################################
# score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  n_eff=ncvar_get(nc,"N_eff")
  basins=ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  n_eff.negis=c()
  for(t in 1:length(time)){
   n_eff.negis[t]=mean(n_eff[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*n_eff.negis
  print(i)
}

mean.n_eff.negis=c()
for(t in 1:length(time)){
  mean.n_eff.negis[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-neff-negis.txt",sep="")
write.table(mean.n_eff.negis, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  n_eff=ncvar_get(nc,"N_eff")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  for(t in 1:length(time)){
    n_eff.negis[t]=mean(n_eff[,,t][basins>9.0 & basins<9.3])
  }
  pd=length(time)
  nc_close(nc)
  nom[i,]=ss.norm[i]*(n_eff.negis-mean.n_eff.negis)^2
  print(i)
}

sd.n_eff.negis=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.n_eff.negis[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-neff-negis.txt",sep="")
write.table(sd.n_eff.negis, file.out, quote = F, sep = " ", row.names=F)


