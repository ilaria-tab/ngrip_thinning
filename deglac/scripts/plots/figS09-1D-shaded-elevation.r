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

################################################################
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


######################################################################################
# s, Hice, z_bed

plot.out=paste(work.fldr,"/Figures/figS09-1D-shaded-elevation.png", sep="")
png(plot.out, width=5, height=6, units="in", res=100)
#dev.new(width=5, height=6, units="in", res=100)
#par(mfrow=c(6,1), mar=c(0,4,0,4), oma=c(4,3,1,4))
layout(matrix(c(1,2,3),nrow=3,byrow=T), heights=c(1,1,1))
#layout.show(n=6)
par(mar=c(0,4,0,4), oma=c(4,3,1,2))

# S
at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
#at.y=seq(-40,240, by=40)
at.y=seq(-50,150, by=50)
xlim=c(-15,0)
#ylim=c(-40,240)
ylim=c(-50,150)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#E2F3EB",50))
#polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#D4EDE0",80))
#lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#9AD5B8",30), lty=3,lwd=2)
#lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#9AD5B8",30), lty=3, lwd=2)
#lines(t.ds.v09/1000, ds.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#9AD5B8",50), lty=1, lwd=2)
#polygon(c(ty, rev(ty)), c(mean.dzs-2*sd.dzs, rev(mean.dzs+2*sd.dzs)), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c(mean.dzs-1*sd.dzs, rev(mean.dzs+1*sd.dzs)), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.dzs, col="#2200B8", lwd=3)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NGRIP S (m)",side=2,line=3.5,cex=1)
text(-15,3150,"a",font=2, cex=2)

# z_bed
at.y=seq(-75,25, by=25)
ylim=c(-75,25)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-2*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+2*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c((mean.zb.ngrip-mean.zb.ngrip[nbands])-1*(sd.zb.ngrip), rev((mean.zb.ngrip-mean.zb.ngrip[nbands])+1*(sd.zb.ngrip))), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.zb.ngrip-mean.zb.ngrip[nbands], lwd=3, col="#2200B8")
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NGRIP bedrock uplift (m)",side=4,line=3.5, cex=1)

# H_ice
at.y=seq(-100,200, by=50)
ylim=c(-100,200)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#grid()
#polygon(c(ty, rev(ty)), c((mean.h.ngrip-mean.h.ngrip[nbands])-2*(sd.h.ngrip), rev((mean.h.ngrip-mean.h.ngrip[nbands])+2*(sd.h.ngrip))), border = NA, col=alpha.col("#305AF2",30))
polygon(c(ty, rev(ty)), c((mean.h.ngrip-mean.h.ngrip[nbands])-1*(sd.h.ngrip), rev((mean.h.ngrip-mean.h.ngrip[nbands])+1*(sd.h.ngrip))), border = NA, col=alpha.col("#305AF2",50))
lines(ty, mean.h.ngrip-mean.h.ngrip[nbands], lwd=3, col="#2200B8")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("NGRIP Ice thickness (m)",side=2,line=3.5, cex=1)

axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr)",side=1,line=2.75,cex=1)


dev.off()


