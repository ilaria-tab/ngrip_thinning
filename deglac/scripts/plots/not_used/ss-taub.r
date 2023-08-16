# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,9],breaks = nr.best))]

# col alpha
alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}


###############################################################################################################
# single scores vs taub

plot.out=paste(work.fldr,"/Figures/sh-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)

#S_H
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S H_PD")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.h[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}

dev.off()

# S_U
plot.out=paste(work.fldr,"/Figures/su-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)

plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S U_PD")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.u[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S_ area_pd
plot.out=paste(work.fldr,"/Figures/sarea-pd-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)

plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S AREA_PD")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.area.pd[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S_NGRIP_PD
plot.out=paste(work.fldr,"/Figures/s-ngrip-pd-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S NGRIP_PD")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.ngrip[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S_NGRIP_DEG
plot.out=paste(work.fldr,"/Figures/s-ngrip-deg-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S NGRIP_DEG")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.cores[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S_79_GL
plot.out=paste(work.fldr,"/Figures/s-79-gl-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S 79_GL")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.79.gl[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S H_NE
plot.out=paste(work.fldr,"/Figures/s-h-ne-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S H_NE")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.h.ne[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S U_NE
plot.out=paste(work.fldr,"/Figures/s-u-ne-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S U_NE")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.u.ne[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

# S_AREA_NE
plot.out=paste(work.fldr,"/Figures/s-area-ne-taub.png", sep="")
png(plot.out, width = 4, height = 4, units = "in", res=100)
#dev.new(width = 10, height = 10, units="in")
xlim=c(50000,70000)
ylim=c(0,1)
plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="S AREA_NE")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, s.area.ne[i], pch=21,cex=1.3, col="black", bg=mycol[i])
}

dev.off()

