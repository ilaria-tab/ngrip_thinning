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
mean.smb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/mean-smb-ngrip.txt",sep=""), skip=1)
mean.smb.ngrip=unlist(mean.smb.ngrip, use.names=FALSE)
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
sd.smb.ngrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-smb-ngrip.txt",sep=""), skip=1)
sd.smb.ngrip=unlist(sd.smb.ngrip, use.names=FALSE)

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
sd.beta.egrip=read.table(paste(work.fldr,"/scripts/score-weighting/sd-beta-egrip.txt",sep=""), skip=1)
sd.beta.egrip=unlist(sd.beta.egrip, use.names=FALSE)

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

plot.out=paste(work.fldr,"/Figures/figS10-lag-thinning-tau.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in", res=100)
layout(matrix(c(1,2,3),nrow=3,byrow=T), heights=c(1,1,1))
#layout.show(n=6)
par(mar=c(0,4,0,4), oma=c(4,3,1,2))

# BETA at NEGIS
at.x=seq(-15,-6,by=2.5)
label.x=seq(15,6,by=-1)
#at.y=seq(-40,240, by=40)
#at.y=seq(30000,70000, by=10000)
at.y=seq(0,8e5, by=1e5)
#label.y=seq(30,70, by=10)
label.y=seq(0,8e5, by=1e5)
xlim=c(-15,-6)
#ylim=c(-40,240)
ylim=c(0,8e5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
#polygon(c(ty, rev(ty)), c(mean.taub.egrip-2*sd.taub.egrip, rev(mean.taub.egrip+1*sd.taub.egrip)), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c(mean.beta.egrip-1*sd.beta.egrip, rev(mean.beta.egrip+1*sd.beta.egrip)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.beta.negis, lwd=5, col="#2200B8")
axis(2,at=at.y,labels=label.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NEGIS Basal fric coeff (kPa)",side=2,line=4, cex=1)
#abline(v=-12, col="black")
#rect(-11.9,28000,-11.7,70000,col=alpha.col("magenta",30), border=NA)
abline(v=-9.8, col=alpha.col("yellow",50), lty=2, lwd=3)
abline(v=-11.4, col=alpha.col("magenta",50), lty=2, lwd=3)
abline(v=-14, col=alpha.col("orange", 50), lwd=3)

# BETA in EGRIP
at.x=seq(-15,-6,by=2.5)
label.x=seq(15,6,by=-1)
#at.y=seq(-40,240, by=40)
#at.y=seq(30000,70000, by=10000)
at.y=seq(6e5,1.4e6, by=1e5)
#label.y=seq(30,70, by=10)
label.y=seq(6e5,1.4e6, by=1e5)
xlim=c(-15,-6)
#ylim=c(-40,240)
ylim=c(6e5,1.4e6)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
#polygon(c(ty, rev(ty)), c(mean.taub.egrip-2*sd.taub.egrip, rev(mean.taub.egrip+1*sd.taub.egrip)), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c(mean.beta.egrip-1*sd.beta.egrip, rev(mean.beta.egrip+1*sd.beta.egrip)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.beta.egrip, lwd=5, col="#2200B8")
axis(2,at=at.y,labels=label.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("EGRIP Basal fric coeff (kPa)",side=2,line=4, cex=1)
#abline(v=-12, col="black")
#rect(-11.9,28000,-11.7,70000,col=alpha.col("magenta",30), border=NA)
abline(v=-9.8, col=alpha.col("yellow",50), lwd=3, lty=2)
abline(v=-11.4, col=alpha.col("magenta",50), lwd=3)
abline(v=-14, col=alpha.col("orange",50), lwd=3, lty=2)

## UB in NGRIP
#at.x=seq(-15,-6,by=1)
#label.x=seq(15,6,by=-1)
##at.y=seq(-40,240, by=40)
##at.y=seq(15000,35000, by=10000)
#at.y=seq(0,0.1,by=0.02)
##label.y=seq(15,35, by=10)
#xlim=c(-15,-6)
##ylim=c(-40,240)
##ylim=c(15000,35000)
#ylim=c(0,0.1)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
##polygon(c(ty, rev(ty)), c(mean.taub.ngrip-2*sd.taub.ngrip, rev(mean.taub.ngrip+1*sd.taub.ngrip)), border = NA, col=alpha.col("#305AF2",30))
##polygon(c(ty, rev(ty)), c(mean.ub.ngrip-1*sd.ub.ngrip, rev(mean.ub.ngrip+1*sd.ub.ngrip)), border = NA, col=alpha.col("#305AF2",50))
#lines(ty, mean.ub.ngrip, lwd=3, col="#2200B8")
#axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("NGRIP Ub (m/yr)",side=4,line=3.5, cex=1)
##abline(v=-12, col="black")
##abline(v=-11, col="magenta")
##rect(-12.4,15000,-11.6,35000,col=alpha.col("magenta",30), border=NA)
#rect(-11.9,14000,-11.7,38000,col=alpha.col("magenta",30), border=NA)
#rect(-10.7,14000,-10.5,38000,col=alpha.col("orange",30), border=NA)


# HICE in NGRIP
at.x=seq(-15,-6,by=1)
label.x=seq(15,6,by=-1)
#at.y=seq(-40,240, by=40)
at.y=seq(3100,3300, by=100)
xlim=c(-15,-6)
#ylim=c(-40,240)
ylim=c(3100,3300)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
#polygon(c(ty, rev(ty)), c(mean.h.ngrip-2*sd.h.ngrip, rev(mean.h.ngrip+1*sd.h.ngrip)), border = NA, col=alpha.col("#305AF2",30))
#polygon(c(ty, rev(ty)), c(mean.h.ngrip-1*sd.h.ngrip, rev(mean.h.ngrip+1*sd.h.ngrip)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.h.ngrip, lwd=5, col="#2200B8")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("H NGRIP (m)",side=2,line=4, cex=1)
#rect(-12.4,3000,-11.6,3200,col=alpha.col("magenta",30), border=NA)
#rect(-11.4,3000,-10.6,3200,col=alpha.col("orange",30), border=NA)
#rect(-11.9,2900,-11.7,3300,col=alpha.col("magenta",30), border=NA)
#rect(-10.7,2900,-10.5,3300,col=alpha.col("orange",30), border=NA)
#rect(-9.8,2900,-9.5,3300,col=alpha.col("yellow",30), border=NA)
abline(v=-9.8, col=alpha.col("yellow",50), lwd=3)
abline(v=-11.4, col=alpha.col("magenta",50), lwd=3, lty=2)
abline(v=-14, col=alpha.col("orange",50), lwd=3, lty=2)

axis(1,at=at.x, labels=label.x, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr ago)",side=1,line=2.75,cex=1)

dev.off()


