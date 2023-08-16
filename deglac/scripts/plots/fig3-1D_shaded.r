library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

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
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

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
#mean.smb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-smb-ngrip.txt",sep=""), skip=1)
#mean.smb.ngrip=unlist(mean.smb.ngrip, use.names=FALSE)
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
#sd.smb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-smb-ngrip.txt",sep=""), skip=1)
#sd.smb.ngrip=unlist(sd.smb.ngrip, use.names=FALSE)

mean.ub.negis=read.table(paste(work.fldr,"/scripts/score-weighting/mean-ub-negis.txt",sep=""), skip=1)
mean.ub.negis=unlist(mean.ub.negis, use.names=FALSE)
mean.beta.negis=read.table(paste(work.fldr,"/scripts/score-weighting/mean-beta-negis.txt",sep=""), skip=1)
mean.beta.negis=unlist(mean.beta.negis, use.names=FALSE)
mean.ub.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-ub-egrip.txt",sep=""), skip=1)
mean.ub.egrip=unlist(mean.ub.egrip, use.names=FALSE)
mean.beta.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-beta-egrip.txt",sep=""), skip=1)
mean.beta.egrip=unlist(mean.beta.egrip, use.names=FALSE)
mean.ub.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-ub-ngrip.txt",sep=""), skip=1)
mean.ub.ngrip=unlist(mean.ub.ngrip, use.names=FALSE)
mean.um.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-um-ngrip.txt",sep=""), skip=1)
mean.um.ngrip=unlist(mean.um.ngrip, use.names=FALSE)

mean.taub.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-taub-ngrip.txt",sep=""), skip=1)
mean.taub.ngrip=unlist(mean.taub.ngrip, use.names=FALSE)
sd.taub.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-taub-ngrip.txt",sep=""), skip=1)
sd.taub.ngrip=unlist(sd.taub.ngrip, use.names=FALSE)
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
#ty=ty-1.95

plot.out=paste(work.fldr,"/Figures/fig3-1D-shaded.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in", res=100)
layout(matrix(c(1,2,3,4,5,6),nrow=6,byrow=T), heights=c(2,1.2,1.2,1.5,1.5,1.2))
par(mar=c(0,4,0,4), oma=c(4,3,1,2))

# S
at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
at.y=seq(2850,3150, by=50)
xlim=c(-15,0)
ylim=c(2850,3150)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)), border = NA, col=alpha.col("#E2F3EB",50))
#polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#D4EDE0",80))
#lines(t.ds.l13/1000, s.ngrip1.l13, col=alpha.col("#9AD5B8",30), lty=3,lwd=2)
#lines(t.ds.l13/1000, s.ngrip2.l13, col=alpha.col("#9AD5B8",30), lty=3, lwd=2)
#lines(t.ds.v09/1000, s.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#9AD5B8",50), lty=1, lwd=2)

#rect(-10.7,2800,-9.7,3150, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("black",30), lty=1, lwd=1)
#abline(v=-9.7, col=alpha.col("black",30), lty=1, lwd=1)
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)), border = NA, col=alpha.col("#C5E7D4",50))
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#B7E1C9",80))
#polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)), border = NA, col=alpha.col("#C7F9CC",50))
#polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#80ED99",80))
lines(t.ds.l13/1000, s.ngrip1.l13, col=alpha.col("#7EC99E",50), lty=3,lwd=3)
lines(t.ds.l13/1000, s.ngrip2.l13, col=alpha.col("#7EC99E",50), lty=3, lwd=3)
lines(t.ds.v09/1000, s.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#57CC99",80), lty=1, lwd=4)

#polygon(c(ty, rev(ty)), c(mean.dzs-2*sd.dzs+2920, rev(mean.dzs+2*sd.dzs+2920)), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c(mean.dzs-1*sd.dzs+2920, rev(mean.dzs+1*sd.dzs+2920)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.dzs+2920, col="#2200B8", lwd=4)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NGRIP S (m)",side=2,line=4.5,cex=1)
text(-15,3150,"a",font=2, cex=2)

#abline(v=-10.6, col=alpha.col("grey",40), lwd=7)
#abline(v=-9.8,col=alpha.col("grey",40), lwd=7 )

# U
at.y=seq(0,6, by=2)
ylim=c(0,6)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#polygon(c(ty, rev(ty)), c(log10(mean.us.negis-1*sd.us.negis), rev(log10(mean.us.negis+1*sd.us.negis))), border = NA, col=alpha.col("#24B34C",30))
#lines(ty, mean.us.negis, col="#24B34C", lwd=3, lty=1)
#polygon(c(ty, rev(ty)), c(mean.us.egrip-mean.us.egrip[length(ty)]-1*sd.us.egrip, rev(mean.us.egrip-mean.us.egrip[length(ty)]+1*sd.us.egrip)), border = NA, col=alpha.col("#24B34C",20))
#lines(ty, mean.us.egrip-mean.us.egrip[length(ty)], col="#24B34C", lwd=3, lty=1)
#grid()
#rect(-10.7,-2,-9.7,8.5, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("green",80), lty=3, lwd=2)
#abline(v=-9.7, col=alpha.col("steelblue",80), lty=3, lwd=2)
#polygon(c(ty, rev(ty)), c(mean.us.ngrip-2*sd.us.ngrip, rev(mean.us.ngrip+2*sd.us.ngrip)), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c(mean.us.ngrip-1*sd.us.ngrip, rev(mean.us.ngrip+1*sd.us.ngrip)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.us.ngrip, col="#2200B8", lwd=4, lty=1)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext(expression("NGRIP U"["s"]*" (m/yr)"),side=4,line=3.2,cex=1)
text(0,5.5,"b",font=2, cex=2)
#abline(v=-10.7, col=alpha.col("grey",40), lwd=7)

# Taub
#at.y=seq(20000,100000, by=20000)
#ylim=c(20000,100000)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#polygon(c(ty, rev(ty)), c(mean.taub.negis-1*sd.taub.negis, rev(mean.taub.negis+1*sd.taub.negis)), border = NA, col=alpha.col("#24B34C",20))
#lines(ty, mean.taub.negis, lwd=3, col="#24B34C")
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NEGIS Taub (Pa)",side=2,line=5, cex=1)

# z_bed
#at.y=seq(-100,50, by=25)
#ylim=c(-100,50)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-2*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+2*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-1*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+1*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.zb.ngrip-mean.zb.ngrip[nbands], lwd=3, col="#305AF2")
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NGRIP bedrock uplift (m)",side=2,line=5, cex=1)

# Neff
at.y=seq(2e6,3.2e7, by=1e7)
ylim=c(2e6,3.2e7)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#rect(-10.7,1e5,-9.7,3.5e7, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("green",80), lty=3, lwd=2)
#abline(v=-9.7, col=alpha.col("steelblue",80), lty=3, lwd=2)
#polygon(c(ty, rev(ty)), c(mean.n_eff.negis-2*sd.n_eff.negis, rev(mean.n_eff.negis+2*sd.n_eff.negis)), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c(mean.n_eff.negis-1*sd.n_eff.negis, rev(mean.n_eff.negis+1*sd.n_eff.negis)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.n_eff.negis, lwd=4, col="#2200B8")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NEGIS Neff (Pa)",side=2,line=5.5, cex=1)
text(-15,3.0e7,"c",font=2, cex=2)

# 79N margin retreat
at.y=seq(-150,350, by=50)
ylim=c(-150,350)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#rect(-10.7,-180,-9.7,380, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("green",80), lty=3, lwd=2)
#abline(v=-9.7, col=alpha.col("steelblue",80), lty=3, lwd=2)
#polygon(c(ty, rev(ty)), c(mean.gl-2*sd.gl, rev(mean.gl+2*sd.gl)), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c(mean.gl-1*sd.gl, rev(mean.gl+1*sd.gl)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.gl, col="#2200B8", lwd=4)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Margin distance (km)",side=4,line=4, cex=1)
# Larsen curve averaged for all NE border
tab=read.table("/home/titan/gwgk/gwgk005h/work/ngrip/data/Larsen2018_margin_NEGIS_Fig2a-3a.dat",skip=0)
t.lar=tab$V1/1000
gl.lar=tab$V2
lines(t.lar, gl.lar, lty=2, lwd=5, col="#fe7f2d")
abline(h=0, col="grey")
# Data from Bennike Weidick 2001 (Blaso) but taken from Larsen2018 SM
t.min=c(-7441,-7252,-7156,-7151,-7150,-6970,-6930,-6393,-5928,-5941,-5566,-4820,-4790)
t.max=c(-7163,-6980,-6795,-6739,-6736,-6679,-6666,-6020,-5664,-5643,-5262,-4525,-4480)
t.min=t.min/1000
t.max=t.max/1000
#for(i in 1:length(t.max)){
#    arrows(x0=t.min[i],y0=-130,x1=t.max[i],y1=-130, code=3, angle=90, length=0.05, col="grey30", lwd=1)
#}
text(0,330,"d",font=2, cex=2)

# dTann 
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

at.y=seq(-15,7.5, by=2.5)
#xlim=c(-15,0)
ylim=c(-15,7.5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#rect(-10.7,-18,-9.7,10, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("green",80), lty=3, lwd=2)
#abline(v=-9.7, col=alpha.col("steelblue",80), lty=3, lwd=2)
polygon(c(t.tann.l17/1000, rev(t.tann.l17/1000)), c(tann.agas.max.l17, rev(tann.agas.min.l17)), border = NA, col=alpha.col("#7400B8",30))
lines(t.tann.l17/1000,tann.agas.l17, lwd=2, col="#7400B8")
lines(t.b18/1000, dtann.b18.79, lwd=2, col="#F5007E")
#lines(t.b18/1000, dtann.b18.ngrip, lwd=1.5, col="#BA208E")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Tann anomaly (°C)",side=2,line=4.5,cex=1)
text(-15,7,"e",font=2, cex=2)


# Pr fraction
at.y=seq(0.2,1.2,by=0.2)
ylim=c(0.2,1.2)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#rect(-10.7,0.2,-9.7,1.4, border = NA, col=alpha.col("grey",60))
#abline(v=-10.7, col=alpha.col("green",80), lty=3, lwd=2)
#abline(v=-9.7, col=alpha.col("steelblue",80), lty=3, lwd=2)
lines(t.b18/1000, prann.b20.ngrip/prann.b20.ngrip[length(t.b18)], col="#342E37", lwd=3)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NGRIP Prec. fraction",side=4,line=4, cex=1)
#abline(v=-9.3, col="pink", lwd=2)
text(0,1.12,"f",font=2, cex=2)

axis(1,at=at.x, labels=label.x, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr ago)",side=1,line=2.75,cex=1)


dev.off()


######################################################################################
## s, Hice, z_bed
#
#plot.out=paste(work.fldr,"/Figures/1D-shaded-elevation.png", sep="")
#png(plot.out, width=5, height=6, units="in", res=100)
##dev.new(width=5, height=6, units="in", res=100)
##par(mfrow=c(6,1), mar=c(0,4,0,4), oma=c(4,3,1,4))
#layout(matrix(c(1,2,3),nrow=3,byrow=T), heights=c(1,1,1))
##layout.show(n=6)
#par(mar=c(0,4,0,4), oma=c(4,3,1,2))
#
## S
#at.x=seq(-15,0,by=2.5)
#label.x=seq(15,0,by=-2.5)
##at.y=seq(-40,240, by=40)
#at.y=seq(-50,150, by=50)
#xlim=c(-15,0)
##ylim=c(-40,240)
#ylim=c(-50,150)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
##grid()
##polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#E2F3EB",50))
##polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#D4EDE0",80))
##lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#9AD5B8",30), lty=3,lwd=2)
##lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#9AD5B8",30), lty=3, lwd=2)
##lines(t.ds.v09/1000, ds.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#9AD5B8",50), lty=1, lwd=2)
##polygon(c(ty, rev(ty)), c(mean.dzs-2*sd.dzs, rev(mean.dzs+2*sd.dzs)), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c(mean.dzs-1*sd.dzs, rev(mean.dzs+1*sd.dzs)), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.dzs, col="#2200B8", lwd=3)
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NGRIP S (m)",side=2,line=3.5,cex=1)
#text(-15,3150,"a",font=2, cex=2)
#
## z_bed
#at.y=seq(-75,25, by=25)
#ylim=c(-75,25)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
##grid()
##polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-2*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+2*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-1*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+1*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.zb.ngrip-mean.zb.ngrip[nbands], lwd=3, col="#2200B8")
#axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NGRIP bedrock uplift (m)",side=4,line=3.5, cex=1)
#
## H_ice
#at.y=seq(-100,200, by=50)
#ylim=c(-100,200)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
##grid()
##polygon(c(ty, rev(ty)), c((mean.h.ngrip-mean.h.ngrip[nbands])-2*(sd.h.ngrip), rev((mean.h.ngrip-mean.h.ngrip[nbands])+2*(sd.h.ngrip))), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c((mean.h.ngrip-mean.h.ngrip[nbands])-1*(sd.h.ngrip), rev((mean.h.ngrip-mean.h.ngrip[nbands])+1*(sd.h.ngrip))), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.h.ngrip-mean.h.ngrip[nbands], lwd=3, col="#2200B8")
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NGRIP Ice thickness (m)",side=2,line=3.5, cex=1)
#
#axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1.3)
#mtext("Time (kyr)",side=1,line=2.75,cex=1)
#
#
#dev.off()
#



#####################################################################################
##plot.out=paste(work.fldr,"/Figures/map-1D.png", sep="")
##png(plot.out, width=4, height=5, units="in", res=100)
###dev.new(width=4, height=5, units="in")
##par(mar=c(1,1,1,4))
#
##file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
#file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
#zs=raster(file, varname="z_srf")
#zs[zs<0]=NA
#
#file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa-negis.nc"
#basin=raster(file, varname="basin")
##basin[basin!=9]=NA
#negis=basin
##negis[negis!=9]=NA
#
#topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
#htopo=raster(topo,varname="H_ice")
#mask=raster(topo, varname="mask")
#zb=raster(topo,varname="z_bed")
#nc=nc_open(topo)
#lat2D=raster(topo,varname="lat2D")
#lon2D=raster(topo,varname="lon2D")
##x2D=ncvar_get(nc,"x2D")
##y2D=ncvar_get(nc,"y2D")
#xc=ncvar_get(nc,"xc")
#yc=ncvar_get(nc,"yc")
#nc_close(nc)
#
#file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17-paleonegis-1sd.nc"
#nc=nc_open(file)
#zb.negis=raster(file,varname="z_bed")
#d=zb.negis-zb
#d[d>=0]=NA
#
## NE transect-> 79N, Blaso
##xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
##ys=c(-1322,-698) # y coord
#xs=c(200,744)
#ys=c(-1354,-690)
#fit=lm(ys~xs)
#x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
#y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
#dat.new=data.frame(x.line, y.line)
#
#xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698) # y coord
#fit=lm(ys~xs)
#x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
#y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
#dat.old=data.frame(x.line, y.line)
#
#plot.out=paste(work.fldr,"/Figures/map-1D.png", sep="")
#png(plot.out, width=4, height=5, units="in", res=100)
##dev.new(width=4, height=5, units="in")
#par(mar=c(1,1,1,1))
#
#plot(zs, col=c("#dee2e6"), box=F,axes=F,legend=F)
#contour(zs, add=T, drawlabels=F, col="grey60")
##plot(negis, col=alpha.col("grey20",30), add=T,legend=F)
#plot(d, add=T, zlim=c(-150,0), )
#contour(negis, level=c(9),add=T, drawlabels=F,lwd=1,col="grey20")
#contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
#  contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
#points(xc[ngrip.i],yc[ngrip.j],bg="grey20",pch=23, cex=1,lwd=1)
##points(xc[egrip.i],yc[egrip.j],col="#6CC25B", pch=20, cex=2)
##points(xc[149],yc[305],bg="grey40", pch=21, cex=1, lwd=1)
##segments(x0=200, y0=-1354, x1=744 , y1=-690,col="orange", lty=1, lwd =2) # new segment
#segments(x0=288, y0=-1322, x1=640 , y1=-698,col="magenta", lty=1, lwd =2) # old segment
##segments(x0=304, y0=-1322, x1=640 , y1=-698,col="magenta", lty=1, lwd =2) # old segment
##segments(x0=288, y0=-1322, x1=704, y1=-698, col="magenta", lty=1, lwd=2)
##segments(x0=272, y0=-1322, x1=704, y1=-698, col="magenta", lty=1, lwd=2)
#points(xc[149],yc[305],bg="grey40", pch=21, cex=1, lwd=1)
#text(0,-1700,"NGRIP",cex=0.8,font=2)
#text(140,-1400,"NEGIS",cex=0.8,font=2)
#text(650,-1000,"79N",cex=0.8,font=2)
##points(xc[156],yc[306],col="#1098F7", pch=20, cex=2)
#dev.off()


