# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("thinning")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

ss.norm=b.best[,17]

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
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  zs=ncvar_get(nc,"z_srf")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd] 
  nom[i,]=ss.norm[i]*dzs.ngrip
  #Sys.sleep(0.5)
  print(i) 
}

mean.dzs=c()
for(t in 1:length(time)){
  mean.dzs[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-dzs.txt",sep="")
write.table(mean.dzs, file.out, quote = F, sep = " ", row.names=F)

# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  zs=ncvar_get(nc,"z_srf")
  basins =ncvar_get(nc,"basins")
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
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  nom[i,]=ss.norm[i]*(dzs.ngrip-mean.dzs)^2
  #Sys.sleep(2)
  print(i)
}

sd.dzs=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.dzs[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-dzs.txt",sep="")
write.table(sd.dzs, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/scripts/score-weighting/time.txt",sep="")
write.table(time, file.out, quote = F, sep = " ", row.names=F)

#i=1
#file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
#nc = nc_open(file)
#xc=ncvar_get(nc, "xc")
#yc= ncvar_get(nc, "yc")
#zs=ncvar_get(nc,"z_srf")
#basins =ncvar_get(nc,"basins")
#time=ncvar_get(nc,"time")
#time=time/1000
#pd=length(time)
#nc_close(nc)
#zs.ngrip=zs[ngrip.i,ngrip.j,]
#best.dzs=zs.ngrip-zs.ngrip[pd]


#dev.new()
#plot(time, mean.dzs, type="l", xlim=c(-15,0), ylim=c(-40,200))
#lines(time, mean.dzs+sd.dzs, col="grey")
#lines(time, mean.dzs-sd.dzs, col="grey")
#lines(time, best.dzs, col="red")


