library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

# load precip anomaly from B20 (highP)
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18-TAS-B20-highP.nc"
nc=nc_open(file)
dtas.b18=ncvar_get(nc, "tas")
pr.b20=ncvar_get(nc,"pr")
t.b18=ncvar_get(nc, "time")
nc_close(nc)
pr.b20.ngrip=pr.b20[ngrip.i,ngrip.j,,]
prann.b20.ngrip=c()
for(t in 1:length(t.b18)){
  prann.b20.ngrip[t]=mean(pr.b20.ngrip[1:12,t])
}

# load present precipitation (BOX2013)
#file="/home/itabone/data/ice_data/Greenland/GRL-8KM/GRL-16KM_BOX2013_monthly_1850-2000-corr.nc"
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BOX2013_monthly_t2m-1981-2010-pr-1850-2000-corr.nc"
nc=nc_open(file)
pr.box=ncvar_get(nc,"pr")
pr.box.ngrip=pr.box[ngrip.i,ngrip.j,]
prann.box.ngrip=mean(pr.box.ngrip) # mm/d
prann.box.ngrip=prann.box.ngrip * 10^(-3)*365

# pr
prann.ngrip=prann.b20.ngrip*prann.box.ngrip

# load simulations
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

# prec 0.5 simulation
file=paste(out.fldr,"/tests/ngrip-prec/305/yelmo2D.nc",sep="")
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
prann=ncvar_get(nc,"Pr_ann")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip.prec=zs[ngrip.i,ngrip.j,]
dzs.ngrip.prec=zs.ngrip.prec-zs.ngrip.prec[pd]
prann.prec=prann[ngrip.i,ngrip.j,]

########################################################################################
# plot
plot.out=paste(work.fldr,"/Figures/figS07-prec.png", sep="")
png(plot.out, width=5, height=7, units="in", res=100)
#dev.new(width=4, height=7, units="in", res=100)
par(mar=c(3,3,1,1), oma=c(2,1,1,1), mfrow=c(2,1))

at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
at.y=seq(0,0.25, by=0.05)
#at.y=seq(2920,3070, by=50)
xlim=c(-15,0)
ylim=c(0,0.25)

plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
lines(t.g14/1000, acc.ngrip.g14, type="l", col="orange",xlab="Time (kyr)", ylab="NGRIP Accumulation (m w.e./yr)")
lines(t.b18/1000, prann.ngrip, col="black", lwd=2)
lines(t.k14/1000, acc.ngrip.k14, col="magenta", lwd=2)
lines(time, prann.prec, col="black", lwd=4, lty=3)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("NGRIP Accumulation (m w.e./yr)",side=2,line=2.5,cex=1)
axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1)
#mtext("Time (kyr ago)",side=1,line=2.75,cex=1)
legend("bottomright", c("B20*Box2013","Kindler et al., 2014","Gkinis et al., 2014"), bty="n",lty=1, lwd=2, col=c("black","magenta","orange"))

at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
#at.y=seq(2900,3100, by=50)
at.y=seq(0,200, by=50)
#at.y=seq(-100,100,by=50)
xlim=c(-15,0)
#ylim=c(2900,3100)
ylim=c(0,200)
#ylim=c(-100,100)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#C5E7D4",50))
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#B7E1C9",80))
lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#7EC99E",50), lty=3,lwd=3)
lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#7EC99E",50), lty=3, lwd=3)
lines(t.ds.v09/1000, ds.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#57CC99",80), lty=1, lwd=4)
lines(time,dzs.ngrip.cntrl, col="black", lwd=3, lty=1)
lines(time,dzs.ngrip.prec, col="black", lwd=3, lty=3)
#diff=dzs.ngrip.prec-dzs.ngrip.cntrl
#lines(time, diff, col="black", lwd=3, lty=1)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("NGRIP dS (m)",side=2,line=2.5,cex=1)
axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.75,cex=1)

dev.off()
