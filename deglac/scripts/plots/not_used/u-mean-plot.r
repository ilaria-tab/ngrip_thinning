# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)
library(TeachingDemos)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,7]

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")


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


#######################################################################################################################
# Velocity fields

# mask areas outside greenland for grounding line 
topo.regions="/home/itabone/data/ice_data/Greenland/GRL-16KM/GRL-16KM_REGIONS.nc"
region=raster(topo.regions, varname="mask")  # Ellesmere 1.11

topo="/home/itabone/data/ice_data/Greenland/GRL-16KM/GRL-16KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
mask=raster(topo, varname="mask")
#nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")
#x2D=ncvar_get(nc,"x2D")
#y2D=ncvar_get(nc,"y2D")
#xc=ncvar_get(nc,"xc")
#yc=ncvar_get(nc,"yc")
#nc_close(nc)
# mask h over observed ice cover
htopo[htopo==0]=NA
# mask ellesmere island
htopo[region < 1.12 & region > 1.10]=NA

bands=c(11,18,22,26,30,34,38,42)

colors=colorRampPalette(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
                          "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
                          "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
                          "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
                          "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
                          "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
                          "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                          "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
                          "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
                          "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
                          "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
                          "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
                          "#AF0000", "#9F0000", "#8F0000", "#800000"))

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

#### band=11 (20 ka BP)
iband=11
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  bed = raster(file, varname="z_bed", band=iband)
  bed[bed < 0]=NA
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.lgm=u
mean.u.lgm=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.lgm[mean.u.lgm==0.000] = NA

#band=18 (12 ka BP)
iband=18
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.12ka=u
mean.u.12ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.12ka[mean.u.12ka==0.000] = NA

#band=22 (10 ka BP)
iband=22
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.10ka=u
mean.u.10ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.10ka[mean.u.10ka==0.000] = NA

#band=26 (8 ka BP)
iband=26
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.8ka=u
mean.u.8ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.8ka[mean.u.8ka==0.000] = NA

#band=30 (6 ka BP)
iband=30
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.6ka=u
mean.u.6ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.6ka[mean.u.6ka==0.000] = NA

#band=34 (4 ka BP)
iband=34
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.4ka=u
mean.u.4ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.4ka[mean.u.4ka==0.000] = NA

#band=38 (2 ka BP)
iband=38
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.2ka=u
mean.u.2ka=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.2ka[mean.u.2ka==0.000] = NA

#band=42 (0 ka BP)
iband=42
nom.u=stack(u)
for(i in 2:nr.best){
  #file = paste(out.fldr,"/btq.0.976.cbz0.-255.600.cfngs.0.681.nffdlt.0.037/yelmo2D.nc",sep="")
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  #file = paste(out.fldr,"/itmb.-2.633.itmc.-25.960.btq.0.963.cbz0.-322.800.cfngs.0.331.nffdlt.0.034/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=iband)
  u.ss=ss.norm[i]*u
  nom.u=addLayer(nom.u,u.ss)
}
mean.u.pd=u
mean.u.pd=sum(nom.u)/sum(ss.norm,na.rm=T)
mean.u.pd[mean.u.pd==0.000] = NA

##################################################################################################
# plot
#plot.out=paste(work.fldr,"/Figures/U-hol-mean.png", sep="")
#png(plot.out, width = 10, height = 7, units = "in", res=100)
dev.new(width=10, height=7, units="in")
par(mfrow=c(2,4), mar=c(0.5,0.5,1.5,1), oma=c(1,1,1,3))

# 20 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("20 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.lgm), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)  

# 12 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.12ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# 10 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.10ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# 8 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.8ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# 6 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.6ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# 4 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.4ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# 2 ka
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.2ka), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)

# pd
plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, main=paste("12 ka ago",sep=""), cex.main=1.5)
plot(log10(mean.u.pd), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=1,add=T)
points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)


#if(i==26 | i==42){
# gris.range <- c(0.1, 1000)
# labels=c(0.1,1,10,100,1000)
# plot(log10(u), horizontal=F, legend.only=T, legend.shrink = 0.5, legend.width = 3, zlim=c(log10(0.1),log10(1000)),
#     col=colors(100),
#     axis.args=list(at=c(log10(0.1),log10(1),log10(10),log10(100),log10(1000)),
#                    labels=labels,
#                    cex.axis=1.3),
#       #legend.args=list(text='U (m/yr)', side=3, font=2, line=0.5, cex=0.8)
# )
#}

#dev.off()
