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


#################################################################################################

# load ds, ghf, Tice

# control simulation
file=paste(out.fldr,"/305/yelmo2D.nc",sep="")
nc=nc_open(file)
basins=ncvar_get(nc,"basins")
u=ncvar_get(nc,"uxy_bar")
zs=ncvar_get(nc,"z_srf")
ghf=ncvar_get(nc,"Q_geo")
tice=ncvar_get(nc,"T_ice")
fpmp=ncvar_get(nc,"f_pmp")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)

fpmp.ngrip.cntrl=fpmp[ngrip.i,ngrip.j,]
fpmp.egrip.cntrl=fpmp[egrip.i,egrip.j,]
u.ngrip.cntrl=u[ngrip.i,ngrip.j,]
u.egrip.cntrl=u[egrip.i,egrip.j,]
zs.ngrip.cntrl=zs[ngrip.i,ngrip.j,]
dzs.ngrip.cntrl=zs.ngrip.cntrl-zs.ngrip.cntrl[pd]
ghf.cntrl=mean(ghf[basins<9.5 & basins>9.0])
ghf.egrip.cntrl=ghf[egrip.i,egrip.j]
ghf.ngrip.cntrl=ghf[ngrip.i,ngrip.j]
tb.negis.cntrl=tb.egrip.cntrl=tb.ngrip.cntrl=c()
for(t in 1:length(time)){
  tbed=tice[,,1,t]
  tb.egrip.cntrl[t]=tice[egrip.i,egrip.j,1,t]
  tb.ngrip.cntrl[t]=tice[ngrip.i,ngrip.j,1,t]
  tb.negis.cntrl[t]=mean(tbed[basins<9.5 & basins>9.0])
}

# rate of change of velocity (acceleration) # m/yr/200yr 
a.ngrip.cntrl=c()
for(t in 1:length(time)){
   a.ngrip.cntrl[t]=(u.ngrip.cntrl[t]-u.ngrip.cntrl[t-1])/(time[t]-time[t-1])
}
a.ngrip.cntrl=a.ngrip.cntrl/200 # m/yr^2



# GHF from S04
file=paste(out.fldr,"/tests/ghf_s04/305/yelmo2D.nc",sep="")
nc=nc_open(file)
basins=ncvar_get(nc,"basins")
u=ncvar_get(nc,"uxy_bar")
zs=ncvar_get(nc,"z_srf")
ghf=ncvar_get(nc,"Q_geo")
tice=ncvar_get(nc,"T_ice")
fpmp=ncvar_get(nc,"f_pmp")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)

fpmp.egrip.s04=fpmp[egrip.i,egrip.j,]
fpmp.ngrip.s04=fpmp[ngrip.i,ngrip.j,]
u.ngrip.s04=u[ngrip.i,ngrip.j,]
u.egrip.s04=u[egrip.i,egrip.j,]
zs.ngrip.s04=zs[ngrip.i,ngrip.j,]
dzs.ngrip.s04=zs.ngrip.s04-zs.ngrip.s04[pd]
ghf.s04=mean(ghf[basins<9.5 & basins>9.0])
ghf.ngrip.s04=ghf[ngrip.i,ngrip.j]
tb.negis.s04=tb.egrip.s04=tb.ngrip.s04=c()
for(t in 1:length(time)){
  tbed=tice[,,1,t]
  tb.egrip.s04[t]=tice[egrip.i,egrip.j,1,t]
  tb.ngrip.s04[t]=tice[ngrip.i,ngrip.j,1,t]
  tb.negis.s04[t]=mean(tbed[basins<9.5 & basins>9.0])
}
# rate of change of velocity (acceleration) # m/yr/200yr 
a.ngrip.s04=c()
for(t in 1:length(time)){
   a.ngrip.s04[t]=(u.ngrip.s04[t]-u.ngrip.s04[t-1])/(time[t]-time[t-1])
}
a.ngrip.s04=a.ngrip.s04/200 # m/yr^2


# GHF=50 mW/yr
file=paste(out.fldr,"/tests/ghf_const/305/yelmo2D.nc",sep="")
nc=nc_open(file)
basins=ncvar_get(nc,"basins")
u=ncvar_get(nc,"uxy_bar")
zs=ncvar_get(nc,"z_srf")
ghf=ncvar_get(nc,"Q_geo")
tice=ncvar_get(nc,"T_ice")
fpmp=ncvar_get(nc,"f_pmp")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)

fpmp.egrip.50=fpmp[egrip.i,egrip.j,]
fpmp.ngrip.50=fpmp[ngrip.i,ngrip.j,]
u.ngrip.50=u[ngrip.i,ngrip.j,]
u.egrip.50=u[egrip.i,egrip.j,]
zs.ngrip.50=zs[ngrip.i,ngrip.j,]
dzs.ngrip.50=zs.ngrip.50-zs.ngrip.50[pd]
ghf.50=mean(ghf[basins<9.5 & basins>9.0])
ghf.ngrip.50=ghf[ngrip.i,ngrip.j]
tb.negis.50=tb.egrip.50=tb.ngrip.50=c()
for(t in 1:length(time)){
  tbed=tice[,,1,t]
  tb.egrip.50[t]=tice[egrip.i,egrip.j,1,t]
  tb.ngrip.50[t]=tice[ngrip.i,ngrip.j,1,t]
  tb.negis.50[t]=mean(tbed[basins<9.5 & basins>9.0])
}

# rate of change of velocity (acceleration) # m/yr/200yr
a.ngrip.50=c()
for(t in 1:length(time)){
   a.ngrip.50[t]=(u.ngrip.50[t]-u.ngrip.50[t-1])/(time[t]-time[t-1])
}
a.ngrip.50=a.ngrip.50/200 # m/yr^2



###########################################################################


# plot dS, fpmp

plot.out=paste(work.fldr,"/Figures/figS05-ghf.png", sep="")
png(plot.out, width=5, height=6, units="in", res=100)
#dev.new(width=5, height=6, units="in", res=100)
layout(matrix(c(1,2,3),nrow=3,byrow=T), heights=c(1,1,1))
#layout.show(n=6)
par(mar=c(0,4,0,4), oma=c(4,3,1,2))
#color=c("black","blue","grey60")
color=c("#002642","#840032","#E59500")

# dS
at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
at.y=seq(-40,140, by=40)
#at.y=seq(2920,3070, by=50)
xlim=c(-15,0)
ylim=c(-20,140)
#ylim=c(2920,3070)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
lines(time,dzs.ngrip.cntrl, col=color[1], lwd=3)
lines(time,dzs.ngrip.s04, col=color[3], lwd=3)
lines(time,dzs.ngrip.50, col=color[2], lwd=3)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("NGRIP dS (m)",side=2,line=2.5,cex=1)

legend("topright",legend=c("63 mW/m2 (M18)","50 mW/m2 (const)","45 mW/m2 (s04)"),col=c(color[1],color[2],color[3]), lwd=2)

# fpmp
at.y=seq(0,1, by=0.2)
ylim=c(0,1)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
grid()
#lines(time,tb.egrip.cntrl, col=color[1], lwd=2)
#lines(time,tb.egrip.s04, col=color[2], lwd=2)
#lines(time,tb.egrip.50, col=color[3], lwd=2)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Frac of grid point at pmp",side=4,line=3,cex=1)
#abline(h=273, lty=2, lwd=3, col="red")
#lines(time,tb.egrip.cntrl, col=color[1], lwd=2)

lines(time,fpmp.ngrip.cntrl, col=color[1], lwd=3, lty=1)
lines(time,fpmp.ngrip.s04, col=color[3], lwd=3, lty=1)
lines(time,fpmp.ngrip.50, col=color[2], lwd=3, lty=1)

lines(time,fpmp.egrip.cntrl, col=color[1], lwd=2, lty=3)
lines(time,fpmp.egrip.s04, col=color[3], lwd=2, lty=3)
lines(time,fpmp.egrip.50, col=color[2], lwd=2, lty=3)

# Acc
at.y=seq(-0.002,0.004, by=0.001)
ylim=c(-0.002,0.004)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5, lwd=2)
grid()
#lines(time,u.egrip.cntrl, col=color[1], lwd=2, lty=3)
#lines(time,u.egrip.s04, col=color[3], lwd=2, lty=3)
#lines(time,u.egrip.50, col=color[2], lwd=2, lty=3)
abline(h=0, lty=3, lwd=3, col="grey")

lines(time,a.ngrip.cntrl, col=color[1], lwd=3, lty=1)
lines(time,a.ngrip.s04, col=color[3], lwd=3, lty=1)
lines(time,a.ngrip.50, col=color[2], lwd=3, lty=1)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Acceleration (m/a^2)",side=2,line=4,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)

dev.off()
