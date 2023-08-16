library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)
library(TeachingDemos)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test6")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test6")
nbands=42 # how many bands does yelmo2D.nc have?

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

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
##########################################################################################################

# plot dS 
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-16,0,by=2)
label.x=seq(16,0,by=-2)
at.y=seq(-50,250, by=50)
xlim=c(-16,0)
ylim=c(-50,250)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
#grid()
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#7FCD79",30))
lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#7FCD79",70), lty=1,lwd=1.5)
lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#7FCD79",70), lty=1, lwd=1.5)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#81AFEE",30))
lines(t.ds.v09/1000, ds.ngrip.v09, type="l", lty=1,xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#81AFEE",70), lwd=1.5)
#dlith=c("1e22","1e23","1e24","1e25")
mycol=c("red","blue")
#j=1
#for(i in dlith){
#  file = paste(out.fldr,"/Dlth.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
#  nc = nc_open(file)
#  xc=ncvar_get(nc, "xc")
#  yc= ncvar_get(nc, "yc")
#  zs=ncvar_get(nc,"z_srf")
#  basins =ncvar_get(nc,"basins")
#  time=ncvar_get(nc,"time")
#  time=time/1000
#  pd=length(time)
#  nc_close(nc)
#  zs.ngrip=zs[ngrip.i,ngrip.j,]
#  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
#  lines(time,  dzs.ngrip, col=mycol[1], lty=j, lwd=1.5)
#  j=j+1 # type of line
#}
j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  zs=ncvar_get(nc,"z_srf")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  lines(time,  dzs.ngrip, col=mycol[2], lty=j, lwd=1.5)
  j=j+1 # type of line
}
file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
zs=ncvar_get(nc,"z_srf")
basins =ncvar_get(nc,"basins")
time=ncvar_get(nc,"time")
time=time/1000
pd=length(time)
nc_close(nc)
zs.ngrip=zs[ngrip.i,ngrip.j,]
dzs.ngrip=zs.ngrip-zs.ngrip[pd]
lines(time,  dzs.ngrip, col="black", lty=j, lwd=2)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP Surface elevation anomaly (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)

####################################################################################################

# rsl change at NGRIP
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(-100,50, by=50)
xlim=c(-16000,0)
ylim=c(-100,50)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  sl.w=sl[ngrip.i,ngrip.j,]
  zb.w=zb[ngrip.i,ngrip.j,]
  nc_close(nc)
  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, lty=j, lwd=1.5, col="blue")
  j=j+1
}

file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
sl=ncvar_get(nc,"z_sl")
zb=ncvar_get(nc,"z_bed")
sl.w=sl[ngrip.i,ngrip.j,]
zb.w=zb[ngrip.i,ngrip.j,]
nc_close(nc)
rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
lines(time, rsl.y, lty=1, lwd=2, col="black")


rsl.lec=c(-80,-60,-40,-10)
rsl.lec1=c(-40,-30,-20,-5)
time.lec=c(-16000,-12000,-8000,-4000)

#plot(time, rsl.y, xlim=c(-18000,0), type="l")
points(time.lec, rsl.lec, pch=21, col="red")
points(time.lec, rsl.lec1, pch=21, col="black")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP RSL (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)

###########################################################################################################

# bedrock elevation
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(-200,-100, by=20)
xlim=c(-16000,0)
ylim=c(-200,-100)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  zb=ncvar_get(nc,"z_bed")
  zb.w=zb[ngrip.i,ngrip.j,]
  nc_close(nc)
  lines(time, zb.w, lty=j, lwd=1.5, col="blue")
  j=j+1
}

file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
zb=ncvar_get(nc,"z_bed")
zb.w=zb[ngrip.i,ngrip.j,]
nc_close(nc)
lines(time, zb.w, lty=1, lwd=2, col="black")

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP zbed (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)

legend("topleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)
#######################################################################################

# surface elevation S
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(2700,3100, by=50)
xlim=c(-16000,0)
ylim=c(2700,3100)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  zs=ncvar_get(nc,"z_srf")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  pd=length(time)
  nc_close(nc)
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  lines(time,  zs.ngrip, col=mycol[2], lty=j, lwd=1.5)
  j=j+1 # type of line
}
file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
zs=ncvar_get(nc,"z_srf")
basins =ncvar_get(nc,"basins")
time=ncvar_get(nc,"time")
pd=length(time)
nc_close(nc)
zs.ngrip=zs[ngrip.i,ngrip.j,]
dzs.ngrip=zs.ngrip-zs.ngrip[pd]
lines(time,  zs.ngrip, col="black", lty=1, lwd=3)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP Surface elevation(m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (yr ago)",side=1,line=2.5,cex=1)

legend("bottomleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)

#######################################################################################

# H ice
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot 
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(2900,3300, by=50)
xlim=c(-16000,0)
ylim=c(2900,3300)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  h=ncvar_get(nc,"H_ice")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  pd=length(time)
  nc_close(nc)
  h.ngrip=h[ngrip.i,ngrip.j,]
  #dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  lines(time,  h.ngrip, col=mycol[2], lty=j, lwd=1.5)
  j=j+1 # type of line
}
file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
h=ncvar_get(nc,"H_ice")
basins =ncvar_get(nc,"basins")
time=ncvar_get(nc,"time")
pd=length(time)
nc_close(nc)
h.ngrip=h[ngrip.i,ngrip.j,]
#dzs.ngrip=zs.ngrip-zs.ngrip[pd]
lines(time,  h.ngrip, col="black", lty=1, lwd=3)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP Hice (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (yr ago)",side=1,line=2.5,cex=1)

legend("bottomleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)

#######################################################################################

# dzb/dt
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

# Plot 
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(-0.002,0.002, by=0.001)
xlim=c(-16000,0)
ylim=c(-0.002,0.002)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  dzbdt=ncvar_get(nc,"dzbdt")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  pd=length(time)
  nc_close(nc)
  dzbdt.ngrip=dzbdt[ngrip.i,ngrip.j,]
  #dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  lines(time,  dzbdt.ngrip, col=mycol[2], lty=j, lwd=1.5)
  j=j+1 # type of line
}
file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
dzbdt=ncvar_get(nc,"dzbdt")
basins =ncvar_get(nc,"basins")
time=ncvar_get(nc,"time")
pd=length(time)
nc_close(nc)
dzbdt.ngrip=dzbdt[ngrip.i,ngrip.j,]
#dzs.ngrip=zs.ngrip-zs.ngrip[pd]
lines(time,  dzbdt.ngrip, col="black", lty=1, lwd=3)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("NGRIP dzb/dt (m/yr)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (yr ago)",side=1,line=2.5,cex=1)

legend("bottomleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)

########################################################################################

# rsl an 79°Glacier
dev.new(width=6, height=4, units="in")
par(mar = c(0,0,0,0), oma = c(3,3,1,1))
par(mar = c(1,2,1,1)) #Add a space 

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij-8km.r")

# North East Bennike&WEidick 2001 (Blaso NORTH)
loc=findij(79.62,-22.50)
t.min=c(-5290,-6900,-5910,-7010,-6910,-6900,-4560,-6920,-6290,-6890,-6160)
t.max=c(-5000,-6730,-5740,-6870,-6800,-6580,-4430,-6800,-6190,-6760,-5920)
rsl=c(10,1,5,27,24,32,9,0.5,1,1,12)

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-16000,0,by=2000)
label.x=seq(16000,0,by=-2000)
at.y=seq(-100,300, by=100)
xlim=c(-16000,0)
ylim=c(-100,300)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

j=1
tau=c("1000","3000","5000","10000","50000")
for(i in tau){
  file = paste(out.fldr,"/t.",i,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]
  nc_close(nc)
  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, lty=j, lwd=1.5, col="blue")
  j=j+1
}

file = paste(out.fldr,"/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
sl=ncvar_get(nc,"z_sl")
zb=ncvar_get(nc,"z_bed")
sl.w=sl[146,300,]
zb.w=zb[146,300,]
nc_close(nc)
rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
lines(time, rsl.y, lty=1, lwd=2, col="black")

for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col="grey30", lwd=1)
    #arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
  }

#plot(time, rsl.y, xlim=c(-18000,0), type="l")
#points(time.lec, rsl.lec, pch=21, col="red")
#points(time.lec, rsl.lec1, pch=21, col="black")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("Blaso RSL (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)
legend("bottomleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)

#################################################################################
# grounding line retreat NE
xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
ys=c(-1322,-698) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)


# plot distance vs time
#plot.out=paste(work.fldr,"/Figures/gldist-vs-time.png", sep="")
#png(plot.out, width = 7, height = 5, units = "in", res=100)
dev.new(width = 7, height = 5, units="in")
xlim=c(-20000,0)
ylim=c(-200,350)
plot(xlim, ylim, type="n", xlab="Time (yr)", ylab="GL distance (km)")

k=1
for(run in tau){
  file = paste(out.fldr,"/t.",run,".itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")

  nc=nc_open(file)
  mask=ncvar_get(nc, "mask_bed")
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time=ncvar_get(nc,"time")

  # where is the GL at the present?
  for(pt in 1:dim(dat)[1]){
    px=round(dat[pt,1])
    py=round(dat[pt,2])
    # from coord to i,j
    i=which(abs(xc-px)==min(abs(xc-px)))
    j=which(abs(yc-py)==min(abs(yc-py)))
    if(length(i)>1 | length(j)>1){
      i=i[1]
      j=j[1]
    }
    if(mask[i,j,length(time)]==4){
      gl.i=i
      gl.j=j
      break
    }else{
      # nothing
    }
  }
  pd.gl=c(gl.i,gl.j)
#  pd.gl=c(73,150)

  # calc distance from pd.gl
  d=c()
  for(t in 1:length(time)){
    for(pt in 1:dim(dat)[1]){
      px=round(dat[pt,1])
      py=round(dat[pt,2])

      # from coord to i,j
      i=which(abs(xc-px)==min(abs(xc-px)))
      j=which(abs(yc-py)==min(abs(yc-py)))

      if(length(i)>1 | length(j)>1){
        i=i[1]
        j=j[1]
      }

      if(mask[i,j,t]==4){
        gl.i=i
        gl.j=j
        break
      }else{
        # nothing
      }
    }
    d[t]=sqrt((xc[pd.gl[1]]-xc[gl.i])^2 + (yc[pd.gl[2]]-yc[gl.j])^2)
    if(gl.i<pd.gl[1] | gl.j<pd.gl[2]){d[t]=-d[t]}
  }
  # plot
  if(run=="3000"){
    lines(time, d, col="black", lwd=2, lty=1)
  }else{
  lines(time, d, col="blue", lwd=1.5, lty=k)
  }
  k=k+1
}
abline(h=0, col="black")

# Data from Larsen 2018
#abline(h=0, col="grey")
t.data=c(-7900)
t.err=c(2690)
for(i in 1:length(t.data)){
  arrows(x0=(t.data[i]-t.err[i]),y0=-50,x1=(t.data[i]+t.err[i]),y1=-50, code=3, angle=90, length=0.05, col="grey40",lwd=1)
}
points(t.data,-50, col="grey40")

# Data from Bennike Weidick 2001 (Blaso) but taken from Larsen2018 SM
t.min=c(-7441,-7252,-7156,-7151,-7150,-6970,-6930,-6393,-5928,-5941,-5566,-4820,-4790)
t.max=c(-7163,-6980,-6795,-6739,-6736,-6679,-6666,-6020,-5664,-5643,-5262,-4525,-4480)
t.min=t.min
t.max=t.max
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=-150,x1=t.max[i],y1=-150, code=3, angle=90, length=0.05, col="grey20", lwd=1)
    #arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
  }

# Larsen curve averaged for all NE border
tab=read.table("/mnt/lustre/scratch/itabone/ilaria/data/Greenland/Larsen2018_groundline_NEGIS.dat",skip=0)
t.lar=tab$V1
gl.lar=tab$V2
lines(t.lar, gl.lar, lty=2, lwd=4, col="purple")

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("Blaso RSL (m)",side=2,line=3,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=2.5,cex=1)
legend("bottomleft",c("tau=1ka","tau=3ka","tau=5ka","tau=10ka","tau=50ka"), lty=c(1,1,3,4,5),col=c("blue","black","blue","blue","blue"), lwd=1.5)

