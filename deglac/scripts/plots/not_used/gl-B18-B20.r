# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# load B20 results
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best.b20=nrow(b.best)
mean.gl.b20=read.table(paste(work.fldr,"/scripts/score-weighting/mean-gl.txt",sep=""), skip=1)
mean.gl.b20=unlist(mean.gl.b20, use.names=FALSE)
sd.gl.b20=read.table(paste(work.fldr,"/scripts/score-weighting/sd-gl.txt",sep=""), skip=1)
sd.gl.b20=unlist(sd.gl.b20, use.names=FALSE)
#time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
#time=unlist(time, use.names=FALSE)
#ty=time


# load B18 results
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best.b18=nrow(b.best)
mean.gl.b18=read.table(paste(work.fldr,"/scripts/score-weighting/mean-gl.txt",sep=""), skip=1)
mean.gl.b18=unlist(mean.gl.b18, use.names=FALSE)
sd.gl.b18=read.table(paste(work.fldr,"/scripts/score-weighting/sd-gl.txt",sep=""), skip=1)
sd.gl.b18=unlist(sd.gl.b18, use.names=FALSE)
time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time


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


##########################################################################################################################

##########################################################################################################################

# dTann B18
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18-TAS-B20-highP.nc"
nc=nc_open(file)
dtas.b18=ncvar_get(nc, "tas")
t.b18=ncvar_get(nc, "time")
nc_close(nc)
dtas.b18.79=dtas.b18[149,305,,]
dtann.b18.79=c()
for(t in 1:length(t.b18)){
  dtann.b18.79[t]=mean(dtas.b18.79[1:12,t])
}
dtann.b18.79=dtann.b18.79-dtann.b18.79[2181] # shift the anomaly wrt to 1750 CE (t[2181])

# dTann B20
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B20-S3T-highP-monthly.nc"
nc=nc_open(file)
dtas.b20=ncvar_get(nc, "tas")
t.b20=ncvar_get(nc, "time")
nc_close(nc)
dtas.b20.79=dtas.b20[149,305,,]
dtann.b20.79=c()
for(t in 1:length(t.b20)){
  dtann.b20.79[t]=mean(dtas.b20.79[1:12,t])
}
dtann.b20.79=dtann.b20.79-dtann.b20.79[397] # shift the anomaly wrt to 1750 CE (t[397])


# plot
plot.out=paste(work.fldr,"/Figures/gl-B18-B20-S3T.png", sep="")
png(plot.out, width=12, height=7, units="in", res=100)
#dev.new(width=12, height=7, units="in", res=100)
#par(mfrow=c(1,2))
par(fig = c(0,1,0,1))
par(mar=c(5,6,0.5,0.5))

# plot GL retreat
at.y=seq(-150,450, by=100)
ylim=c(-150,450)
xlim=c(-14,0)
at.x=seq(-14,0, by=2)
lab.x=seq(14,0, by=-2)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
polygon(c(ty, rev(ty)), c(mean.gl.b20-2*sd.gl.b20, rev(mean.gl.b20+2*sd.gl.b20)), border = NA, col=alpha.col("#ffbd00",20)) #"#f72585",10))
lines(ty, mean.gl.b20,col="#ffbd00", lwd=5)   #col="#f72585", lwd=5)
polygon(c(ty, rev(ty)), c(mean.gl.b18-2*sd.gl.b18, rev(mean.gl.b18+2*sd.gl.b18)), border = NA, col=alpha.col("#390099",20))#"#0496ff",20))
lines(ty, mean.gl.b18,col="#390099", lwd=5) #col="#0496ff", lwd=5)
axis(1,at=at.x,labels=lab.x,col="black",col.axis="black",las=1, cex.axis=1.75)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.75)
mtext("GL distance (km)",side=2,line=3.75, cex=1.75)
mtext("Time (kyr ago)", side=1, line=3.3,cex=1.75)


# Larsen curve averaged for all NE border
tab=read.table("/home/titan/gwgk/gwgk005h/work/ngrip/data/Larsen2018_margin_NEGIS_Fig2a-3a.dat",skip=0)
t.lar=tab$V1/1000
gl.lar=tab$V2
lines(t.lar, gl.lar, lty=2, lwd=7, col="#ff0054")
abline(h=0, col="grey")

text(-14.3,430, "a", cex=2, font=2)

# plot dTann
par(fig = c(0.6,1, 0.5, 1), new = T)
#par(fig = c(0.02,0.45, 0.02, 0.55), new = T)
at.y=seq(-5,7.5, by=2.5)
xlim=c(-12.5,0)
ylim=c(-5,7.5)
at.x=seq(-12.5,0, by=2.5)
lab.x=seq(12.5,0, by=-2.5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
grid()
polygon(c(t.tann.l17/1000, rev(t.tann.l17/1000)), c(tann.agas.max.l17, rev(tann.agas.min.l17)), border = NA, col=alpha.col("grey60",30))
lines(t.tann.l17/1000,tann.agas.l17, lwd=2, col="grey60")
lines(t.b18/1000, dtann.b18.79, lwd=1.5, col="#390099") #"#0496ff")
lines(t.b20/1000, dtann.b20.79, lwd=2, col="#ffbd00") #"#f72585")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.2)
mtext("Tann anomaly (°C)",side=2,line=2.75,cex=1.2)
axis(1,at=at.x,labels=lab.x,col="black",col.axis="black",las=1, cex.axis=1.2)
mtext("Time (kyr ago)", side=1, line=2, cex=1.2)
text(-12.5,7, "b", cex=2, font=2)

dev.off()

