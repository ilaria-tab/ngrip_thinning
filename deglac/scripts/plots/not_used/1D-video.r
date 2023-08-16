library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test9_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test9_deglac")
nbands=87 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,13],breaks = nr.best))]

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
# load mean+ sd values
mean.dzs=read.table(paste(work.fldr,"/scripts/score-weighting/mean-dzs.txt",sep=""), skip=1)
mean.dzs=unlist(mean.dzs, use.names=FALSE)
sd.dzs=read.table(paste(work.fldr,"/scripts/score-weighting/sd-dzs.txt",sep=""), skip=1)
sd.dzs=unlist(sd.dzs, use.names=FALSE)
mean.us.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-us-egrip.txt",sep=""), skip=1)
mean.us.egrip=unlist(mean.us.egrip, use.names=FALSE)
mean.us.negis=read.table(paste(work.fldr,"/scripts/score-weighting/mean-us-negis.txt",sep=""), skip=1)
mean.us.negis=unlist(mean.us.negis, use.names=FALSE)
mean.us.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-us-ngrip.txt",sep=""), skip=1)
mean.us.ngrip=unlist(mean.us.ngrip, use.names=FALSE)
mean.zb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-zb-ngrip.txt",sep=""), skip=1)
mean.zb.ngrip=unlist(mean.zb.ngrip, use.names=FALSE)
mean.h.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-h-ngrip.txt",sep=""), skip=1)
mean.h.ngrip=unlist(mean.h.ngrip, use.names=FALSE)
sd.us.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-us-ngrip.txt",sep=""), skip=1)
sd.us.ngrip=unlist(sd.us.ngrip, use.names=FALSE)
sd.us.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-us-egrip.txt",sep=""), skip=1)
sd.us.egrip=unlist(sd.us.egrip, use.names=FALSE)
sd.us.negis=read.table(paste(work.fldr,"/scripts/score-weighting/sd-us-negis.txt",sep=""), skip=1)
sd.us.negis=unlist(sd.us.negis, use.names=FALSE)
sd.zb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-zb-ngrip.txt",sep=""), skip=1)
sd.zb.ngrip=unlist(sd.zb.ngrip, use.names=FALSE)
sd.h.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-h-ngrip.txt",sep=""), skip=1)
sd.h.ngrip=unlist(sd.h.ngrip, use.names=FALSE)

mean.taub.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-taub-egrip.txt",sep=""), skip=1)
mean.taub.egrip=unlist(mean.taub.egrip, use.names=FALSE)
mean.taub.negis=read.table(paste(work.fldr,"/scripts/score-weighting/mean-taub-negis.txt",sep=""), skip=1)
mean.taub.negis=unlist(mean.taub.negis, use.names=FALSE)
sd.taub.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-taub-egrip.txt",sep=""), skip=1)
sd.taub.egrip=unlist(sd.taub.egrip, use.names=FALSE)
sd.taub.negis=read.table(paste(work.fldr,"/scripts/score-weighting/sd-taub-negis.txt",sep=""), skip=1)
sd.taub.negis=unlist(sd.taub.negis, use.names=FALSE)
mean.n_eff.negis=read.table(paste(work.fldr,"/scripts/score-weighting/mean-neff-negis.txt",sep=""), skip=1)
mean.n_eff.negis=unlist(mean.n_eff.negis, use.names=FALSE)
sd.n_eff.negis=read.table(paste(work.fldr,"/scripts/score-weighting/sd-neff-negis.txt",sep=""), skip=1)
sd.n_eff.negis=unlist(sd.n_eff.negis, use.names=FALSE)
mean.gl=read.table(paste(work.fldr,"/scripts/score-weighting/mean-gl.txt",sep=""), skip=1)
mean.gl=unlist(mean.gl, use.names=FALSE)
sd.gl=read.table(paste(work.fldr,"/scripts/score-weighting/sd-gl.txt",sep=""), skip=1)
sd.gl=unlist(sd.gl, use.names=FALSE)

time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time

# climate forcing 
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18-TAS-B20-highP.nc"
nc=nc_open(file)
dtas.b18=ncvar_get(nc, "tas")
pr.b20=ncvar_get(nc,"pr")
t.b18=ncvar_get(nc, "time")
nc_close(nc)
#dtas.b18.79=dtas.b18[156,306,,]
dtas.b18.79=dtas.b18[149,305,,]
dtas.b18.ngrip=dtas.b18[ngrip.i,ngrip.j,,]
pr.b20.ngrip=pr.b20[ngrip.i,ngrip.j,,]
#dtas.b18.79=dtas.b18[gl.i,gl.j,,]
dtann.b18.79=c()
for(t in 1:length(t.b18)){
  dtann.b18.79[t]=mean(dtas.b18.79[1:12,t])
}
# shift the anomaly wrt to 1750 CE (t[2181])
dtann.b18.79=dtann.b18.79-dtann.b18.79[2181]
dtann.b18.ngrip=c()
for(t in 1:length(t.b18)){
  dtann.b18.ngrip[t]=mean(dtas.b18.ngrip[1:12,t])
}
dtann.b18.ngrip=dtann.b18.ngrip-dtann.b18.ngrip[2181]

prann.b20.ngrip=c()
for(t in 1:length(t.b18)){
  prann.b20.ngrip[t]=mean(pr.b20.ngrip[1:12,t])
}
## shift the anomaly wrt to 1750 CE (t[2181])
#prann.b20.ngrip=prann.b20.ngrip-prann.b20.ngrip[2181]



#### plot

# from 15kyr to 0. One plot every 0.2 kyr = 75 plots. 
# step for Buizert 2018 temperature data (to plot it every 0.2 kyr)
jtemp=20
# step for yelmo 2D data (to plot it every 0.2 kyr)
jy=1
for(i in 1:75){
  plot.out=paste(work.fldr,"/Figures/video/1D/1D-video",i,".png", sep="")
  png(plot.out, width=5, height=10, units="in", res=100)
  #dev.new(width=5, height=10, units="in", res=100)
  #par(mfrow=c(6,1), mar=c(0,4,0,4), oma=c(4,3,1,4))
  layout(matrix(c(1,2,3,4,5),nrow=5,byrow=T), heights=c(1.5,1.5,1.2,1.2,2))
  #layout.show(n=6)
  par(mar=c(0,4,0,4), oma=c(4,3,1,2))

  # dTann
  at.y=seq(-15,7.5, by=2.5)
  xlim=c(-15,0)
  ylim=c(-15,7.5)
  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
  grid()
  dtann.temp=dtann.b18.79[(701):(701+jtemp)]
  t.b18.temp=t.b18[(701):(701+jtemp)]
  #polygon(c(t.tann.l17/1000, rev(t.tann.l17/1000)), c(tann.agas.max.l17, rev(tann.agas.min.l17)), border = NA, col=alpha.col("grey60",30))
  #lines(t.tann.l17/1000,tann.agas.l17, lwd=2, col="grey60")
  lines(t.b18.temp/1000, dtann.temp, lwd=2, col="#F5007E")
  #lines(t.b18/1000, dtann.b18.ngrip, lwd=1.5, col="#BA208E")
  axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
  mtext("Tann anomaly (°C)",side=2,line=4.5,cex=1)
  #text(-15,7,"e",font=2, cex=2)


  # GL retreat
  at.y=seq(-150,350, by=50)
  ylim=c(-150,350)
  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
  #rect(-9.3,-70,-1,-20, col="mistyrose", border=NA)
  grid()
  #polygon(c(ty, rev(ty)), c(mean.gl-2*sd.gl, rev(mean.gl+2*sd.gl)), border = NA, col=alpha.col("#305AF2",30))
  #polygon(c(ty, rev(ty)), c(mean.gl-1*sd.gl, rev(mean.gl+1*sd.gl)), border = NA, col=alpha.col("#305AF2",50))
  ty.temp=ty[(12):(12+jy)]
  mean.gl.temp=mean.gl[(12):(12+jy)]
  sd.gl.temp=sd.gl[(12):(12+jy)]
  polygon(c(ty.temp, rev(ty.temp)), c(mean.gl.temp-2*sd.gl.temp, rev(mean.gl.temp+2*sd.gl.temp)), border = NA, col=alpha.col("#305AF2",30))
  polygon(c(ty.temp, rev(ty.temp)), c(mean.gl.temp-1*sd.gl.temp, rev(mean.gl.temp+1*sd.gl.temp)), border = NA, col=alpha.col("#305AF2",50))
  lines(ty.temp, mean.gl.temp, col="#2200B8", lwd=5)
  axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
  mtext("GL distance (km)",side=4,line=4, cex=1)
  # Larsen curve averaged for all NE border
  tab=read.table("/home/titan/gwgk/gwgk005h/work/ngrip/data/Larsen2018_margin_NEGIS_Fig2a-3a.dat",skip=0)
  t.lar=tab$V1/1000
  gl.lar=tab$V2
  lines(t.lar, gl.lar, lty=2, lwd=3, col="#fe7f2d")
  # Data from Larsen 2018
  abline(h=0, col="grey")

  # Taub
  at.y=seq(20000,120000, by=20000)
  ylim=c(20000,120000)
  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
  grid()
  #polygon(c(ty, rev(ty)), c(mean.taub.negis-1*sd.taub.negis, rev(mean.taub.negis+1*sd.taub.negis)), border = NA, col=alpha.col("#24B34C",20))
  ty.temp=ty[(12):(12+jy)]
  mean.taub.temp=mean.taub.negis[(12):(12+jy)]
  sd.taub.temp=sd.taub.negis[(12):(12+jy)]
  polygon(c(ty.temp, rev(ty.temp)), c(mean.taub.temp-2*sd.taub.temp, rev(mean.taub.temp+2*sd.taub.temp)), border = NA, col=alpha.col("#305AF2",30))
  polygon(c(ty.temp, rev(ty.temp)), c(mean.taub.temp-1*sd.taub.temp, rev(mean.taub.temp+1*sd.taub.temp)), border = NA, col=alpha.col("#305AF2",50))
  lines(ty.temp, mean.taub.temp, lwd=5, col="#2200B8")
  axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
  mtext("NEGIS Taub (Pa)",side=2,line=5, cex=1)

  # U
  at.y=seq(0,18, by=2)
  ylim=c(0,18)
  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
  grid()
  ty.temp=ty[(12):(12+jy)]
  mean.us.temp=mean.us.egrip[(12):(12+jy)]
  sd.us.temp=sd.us.egrip[(12):(12+jy)]
  polygon(c(ty.temp, rev(ty.temp)), c(mean.us.temp-2*sd.us.temp, rev(mean.us.temp+2*sd.us.temp)), border = NA, col=alpha.col("#305AF2",30))
  polygon(c(ty.temp, rev(ty.temp)), c(mean.us.temp-1*sd.us.temp, rev(mean.us.temp+1*sd.us.temp)), border = NA, col=alpha.col("#305AF2",50))
  lines(ty.temp, mean.us.temp, col="#2200B8", lwd=5, lty=1)
  axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
  mtext("EGRIP Us (m/yr)",side=4,line=3.2,cex=1)
  #text(0,7.5,"b",font=2, cex=2)

  # S
  at.x=seq(-15,0,by=2.5)
  label.x=seq(15,0,by=-2.5)
  #at.y=seq(-40,240, by=40)
  at.y=seq(2850,3150, by=50)
  xlim=c(-15,0)
  #ylim=c(-40,240)
  ylim=c(2850,3150)
  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
  grid()
  polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)), border = NA, col=alpha.col("#C5E7D4",50))
  polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#B7E1C9",80))
  lines(t.ds.l13/1000, s.ngrip1.l13, col=alpha.col("#7EC99E",50), lty=3,lwd=2.5)
  lines(t.ds.l13/1000, s.ngrip2.l13, col=alpha.col("#7EC99E",50), lty=3, lwd=2.5)
  lines(t.ds.v09/1000, s.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#7EC99E",80), lty=1, lwd=2.5)

  #polygon(c(ty, rev(ty)), c(mean.dzs-2*sd.dzs+2920, rev(mean.dzs+2*sd.dzs+2920)), border = NA, col=alpha.col("#305AF2",30))
  #polygon(c(ty, rev(ty)), c(mean.dzs-1*sd.dzs+2920, rev(mean.dzs+1*sd.dzs+2920)), border = NA, col=alpha.col("#305AF2",50))
  ty.temp=ty[(12):(12+jy)]
  mean.dzs.temp=mean.dzs[(12):(12+jy)]
  sd.dzs.temp=sd.dzs[(12):(12+jy)]
  polygon(c(ty.temp, rev(ty.temp)), c(mean.dzs.temp-2*sd.dzs.temp+2920, rev(mean.dzs.temp+2*sd.dzs.temp+2920)), border = NA, col=alpha.col("#305AF2",30))
  polygon(c(ty.temp, rev(ty.temp)), c(mean.dzs.temp-1*sd.dzs.temp+2920, rev(mean.dzs.temp+1*sd.dzs.temp+2920)), border = NA, col=alpha.col("#305AF2",50))
  lines(ty.temp, mean.dzs.temp+2920, col="#2200B8", lwd=5)
  axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
  mtext("NGRIP S (m)",side=2,line=4.5,cex=1)
  #text(-15,3150,"a",font=2, cex=2)

  axis(1,at=at.x, labels=label.x, col="black", col.axis="black", las=1, cex.axis=1.3)
  mtext("Time (kyr ago)",side=1,line=2.75,cex=1)

  jy=jy+1
  jtemp=jtemp+20

  dev.off()
}


# Neff
#at.y=seq(2e6,2.3e7, by=1e7)
#ylim=c(2e6,2.3e7)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#polygon(c(ty, rev(ty)), c(mean.n_eff.negis-2*sd.n_eff.negis, rev(mean.n_eff.negis+2*sd.n_eff.negis)), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c(mean.n_eff.negis-1*sd.n_eff.negis, rev(mean.n_eff.negis+1*sd.n_eff.negis)), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.n_eff.negis, lwd=3, col="#2200B8")
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NEGIS Neff (Pa)",side=2,line=5.5, cex=1)
#text(-15,2.2e7,"c",font=2, cex=2)
