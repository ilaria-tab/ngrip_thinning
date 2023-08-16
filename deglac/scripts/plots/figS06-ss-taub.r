# Plots
  
library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,19],breaks = nr.best))]
best.sim=b.best[1,1]


# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

###############################################################################################################
# single scores vs taub

plot.out=paste(work.fldr,"/Figures/figS06-ss-taub.png", sep="")
png(plot.out, width = 15, height = 10, units = "in", res=100)
#dev.new(width = 15, height = 10, units="in")
par(mfrow=c(3,5))
xlim=c(45000,75000)
ylim=c(0,1)

print("Start: area_lgm")
#S_AREA_LGM
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S AREA_LGM")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.area.lgm[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}

file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
nc=nc_open(file)
basins=ncvar_get(nc,"basins")
taub=ncvar_get(nc,"taub")
time=ncvar_get(nc,"time")
nc_close(nc)
  # calc thinning rate
max.i=which(zs.ngrip==max(zs.ngrip))
pd.i=length(time)
min.i = pd.i
#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
taub.hol=taub[,,max.i:min.i]
taub.negis.best=mean(taub.hol[basins>8.5 & basins<9.5])
points(taub.negis.best,s.area.lgm[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)


#S_H_PD
print("Start: h_pd")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S H_PD")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.h[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.h[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)


#S_AREA_PD
print("Start: area_pd")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S AREA_PD")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.area.pd[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.area.pd[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)



#S_U_PD
print("Start: u_pd")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S U_PD")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.u[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.u[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)

#S_NGRIP_PD
print("Start: ngrip_pd")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S NGRIP_PD")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.cores[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.cores[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)




#S_AREA_NE
print("Start: area_ne")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S AREA_NE")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.area.ne[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.area.ne[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)




#S_H_NE
print("Start: h_ne")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S H_NE")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.h.ne[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}

points(taub.negis.best,s.h.ne[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)

#S_U_NE
print("Start: u_ne")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S U_NE")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.u.ne[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.u.ne[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)


#S_NGRIP_DEG
print("Start: ngrip_deg")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S NGRIP_DEG")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.ngrip[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.ngrip[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)



#S_SEDIM_DEG
print("Start: sedim_deg")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S SEDIM_DEG")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.marine[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.marine[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)



#S_OUT_MORAIN_DEG
print("Start: out_morain_deg")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S OUT_MORAIN_DEG")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.outer.coast[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.outer.coast[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)



#S_INN_MORAIN_DEG
print("Start: inn_morain_deg")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S INN_MORAIN_DEG")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.inner.coast[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.inner.coast[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)




#S_79N_DEG
print("Start: 79N_deg")
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S 79N_DEG")
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

#  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])
  points(taub.negis, s.ng[i], pch=21,cex=1.3, col="black", bg=alpha.col(mycol[i],50))
  print(i)
}
points(taub.negis.best,s.ng[best.sim], pch=24, cex=2,col="blue", bg=mycol[best.sim], lwd=3)




dev.off()


