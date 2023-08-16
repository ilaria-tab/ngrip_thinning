library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90  # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs (ngrip thinning)
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,19],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

########################################################################################

# control simulation
file=paste(out.fldr,"/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.cntrl=zs[ngrip.i,ngrip.j,]
dzs.ngrip.cntrl=zs.ngrip.cntrl-zs.ngrip.cntrl[pd]

# case 7 (8kyr)
file=paste(out.fldr,"/tests/switch-case7-8kyr/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.case7=zs[ngrip.i,ngrip.j,]
dzs.ngrip.case7=zs.ngrip.case7-zs.ngrip.case7[pd]

# case 8 (8 kyr)
file=paste(out.fldr,"/tests/switch-case8-8kyr/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.case8=zs[ngrip.i,ngrip.j,]
dzs.ngrip.case8=zs.ngrip.case8-zs.ngrip.case8[pd]

# case 6 (4kyr)
file=paste(out.fldr,"/tests/switch-case6-4kyr/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.case6=zs[ngrip.i,ngrip.j,]
dzs.ngrip.case6=zs.ngrip.case6-zs.ngrip.case6[pd]

# case 7 (4kyr)
file=paste(out.fldr,"/tests/switch-case7-4kyr/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.case7.4kyr=zs[ngrip.i,ngrip.j,]
dzs.ngrip.case7.4kyr=zs.ngrip.case7.4kyr-zs.ngrip.case7.4kyr[pd]

# case 8 (4 kyr)
file=paste(out.fldr,"/tests/switch-case8-4kyr/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.case8.4kyr=zs[ngrip.i,ngrip.j,]
dzs.ngrip.case8.4kyr=zs.ngrip.case8.4kyr-zs.ngrip.case8.4kyr[pd]
###########################################################################

# plot dS

plot.out=paste(work.fldr,"/Figures/figS08-switch.png", sep="")
png(plot.out, width=4, height=10, units="in", res=100)
#dev.new(width=4, height=10, units="in", res=100)
layout(matrix(c(1,2,3,4),nrow=4,byrow=T), heights=c(1,1,1,1))
#layout.show(n=6)
par(mar=c(2,4,1,1), oma=c(4,3,1,1))
#color=c("black","blue","grey60")
color=c("#002642","#840032","#E59500")

at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
at.y=seq(0,140, by=20)
#at.y=seq(2920,3070, by=50)
xlim=c(-15,0)
ylim=c(0,140)

plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
lines(time,dzs.ngrip.cntrl, col=color[1], lwd=3)
lines(time,dzs.ngrip.case7, col=color[3], lwd=3)
lines(time,dzs.ngrip.case8, col=color[2], lwd=3)
lines(time,dzs.ngrip.case6, col=color[1], lwd=3, lty=2)
lines(time,dzs.ngrip.case7.4kyr, col=color[3], lwd=3, lty=2)
lines(time,dzs.ngrip.case8.4kyr, col=color[2], lwd=3, lty=2)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("NGRIP dS (m)",side=2,line=2.5,cex=1)
axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1.3)

legend("topright",legend=c("Cntrl","Case 2","Case 3"),col=c(color[1],color[3],color[2]), lwd=2)
######################################################################

# plot switch case 10 (cntrl)
at.x=c(-15,-8,0)
at.y=c(0,0.3,1)
xlim=c(-15,0)
ylim=c(0,1)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n",yaxt="n", cex.axis=1.5)
segments(-15,1,-8.01,1, col="#01E63A",lwd=6)
segments(-8.01,1,-8,0.3, col="#01E63A", lwd=6)
segments(-8,0.3,0,0, col="#01E63A", lwd=6)
segments(-15,0.3,-8,0.3, col="#ff006e",lwd=6)
segments(-8,0.3,0,1, col="#ff006e", lwd=6)
segments(-15,0.3,0,0.3, col="#ffbe0b", lwd=6, lty=3)
text(-4,0.73,"North",cex=1.2, font=3, srt=43)
text(-4,0.35,"Centre",cex=1.2, font=3, srt=0)
text(-4,0.085,"South",cex=1.2, font=3, srt=-20)
axis(2,at=at.y,labels=c("cf_NEGIS*k","cf_NEGIS","1"),las=3,col="black",col.axis="black",las=1, cex.axis=1.3)
axis(1,at=at.x,labels=c("-15","ts","0"), col="black", col.axis="black", las=1, cex.axis=1.3)
#mtext("Time (kyr ago)",side=1,line=2.75,cex=1.3)
text(-12,0.8,"Cntrl", col=color[1],font=2, cex=1.5)

# plot switch case 7 (4.5 and 8 kyr)
at.x=c(-15,-8,0)
at.y=c(0,0.3,1)
xlim=c(-15,0)
ylim=c(0,1)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n",yaxt="n", cex.axis=1.5)
segments(-15,1,-8.1,1, col="#01E63A",lwd=6)
segments(-8.1,1,-8.05,0.3, col="#01E63A", lwd=6)
segments(-8.05,0.3,0,0, col="#01E63A", lwd=6)
segments(-15,0.3,-7.9,0.3, col="#ff006e",lwd=6)
segments(-7.9,0.3,-7.8,1, col="#ff006e", lwd=6)
segments(-7.8,1,0,1, col="#ff006e", lwd=6)
segments(-15,0.3,0,0.3, col="#ffbe0b", lwd=6, lty=3)
text(-4,0.9,"North",cex=1.2, font=3, srt=0)
text(-4,0.35,"Centre",cex=1.2, font=3, srt=0)
text(-4,0.085,"South",cex=1.2, font=3, srt=-20)
axis(2,at=at.y,labels=c("cf_NEGIS*k","cf_NEGIS","1"),las=3,col="black",col.axis="black",las=1, cex.axis=1.3)
axis(1,at=at.x, labels=c("-15","ts","0"),las=3, col="black", col.axis="black", las=1, cex.axis=1.3)
#mtext("Time (kyr ago)",side=1,line=2.75,cex=1.3)
text(-12,0.8,"Case 2", col=color[3], font=2, cex=1.5)

# plot switch case 8 (4.5 and 8 kyr)
at.x=c(-15,-8,0)
at.y=c(0,0.3,1)
xlim=c(-15,0)
ylim=c(0,1)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n",yaxt="n", cex.axis=1.5)
segments(-15,0.29,-8,0.29, col="#01E63A",lwd=6)
segments(-8,0.29,0,0, col="#01E63A", lwd=6)
segments(-15,0.31,-8,0.31, col="#ff006e",lwd=6)
segments(-8,0.31,0,1, col="#ff006e", lwd=6)
segments(-15,0.3,0,0.3, col="#ffbe0b", lwd=6, lty=3)
text(-4,0.73,"North",cex=1.2, font=3, srt=43)
text(-4,0.35,"Centre",cex=1.2, font=3, srt=0)
text(-4,0.085,"South",cex=1.2, font=3, srt=-20)
axis(2,at=at.y,labels=c("cf_NEGIS*k","cf_NEGIS","1"),las=3,col="black",col.axis="black",las=1, cex.axis=1.3)
axis(1,at=at.x, labels=c("-15","ts","0"),las=3, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr ago)",side=1,line=2.75,cex=1.3)
text(-12,0.8,"Case 3", col=color[2], font=2, cex=1.5)


dev.off()
