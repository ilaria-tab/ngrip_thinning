# Plots

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

##########################################################################################################################
## Plot dS with shaded areas (ensemble mean + sd)
#mean.dzs=read.table(paste(work.fldr,"/scripts/score-weighting/mean-dzs.txt",sep=""), skip=1)
#mean.dzs=unlist(mean.dzs, use.names=FALSE)
#sd.dzs=read.table(paste(work.fldr,"/scripts/score-weighting/sd-dzs.txt",sep=""), skip=1)
#sd.dzs=unlist(sd.dzs, use.names=FALSE)
#time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
#time=unlist(time, use.names=FALSE)
##time=time-1.95
#
#plot.out=paste(work.fldr,"/Figures/dS-NGRIP-shaded.png", sep="")
#png(plot.out, width = 6, height = 4, units = "in", res=100)
##dev.new(width=6, height=4, units="in")
#par(mar = c(0,0,0,0), oma = c(3,3,1,1))
#par(mar = c(1,2,1,1)) #Add a space 
#
## Plot dS
#par(fig=c(0,1,0,1))
#at.x=seq(-16,0,by=2)
#label.x=seq(16,0,by=-2)
#at.y=seq(-50,250, by=50)
#xlim=c(-16,0)
#ylim=c(-50,250)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
##grid()
#polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#7FCD79",30))
#lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#7FCD79",70), lty=1,lwd=1.5)
#lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#7FCD79",70), lty=1, lwd=1.5)
#polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#81AFEE",30))
#lines(t.ds.v09/1000, ds.ngrip.v09, type="l", lty=1,xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#81AFEE",70), lwd=1.5)
#
## shaded ensemble mean + 2*sd
#polygon(c(time, rev(time)), c(mean.dzs-2*sd.dzs, rev(mean.dzs+2*sd.dzs)), border = NA, col=alpha.col("#E95D8C",20))
#lines(time, mean.dzs, col="#D71D5B", lty=1, lwd=2)
##lines(time, mean.dzs, col="#4CB944", lty=1, lwd=2)
#
## shaded ensemble mean + sd
#polygon(c(time, rev(time)), c(mean.dzs-1*sd.dzs, rev(mean.dzs+1*sd.dzs)), border = NA, col=alpha.col("#E95D8C",60))
#lines(time, mean.dzs, col="#D71D5B", lty=1, lwd=2)
##lines(time, mean.dzs, col="#4CB944", lty=1, lwd=2)
#
#file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc",sep="")
#nc = nc_open(file)
#xc=ncvar_get(nc, "xc")
#yc= ncvar_get(nc, "yc")
#zs=ncvar_get(nc,"z_srf")
#basins =ncvar_get(nc,"basins")
#time=ncvar_get(nc,"time")
#time=time/1000
##time=time-1.95
#pd=length(time)
#nc_close(nc)
#zs.ngrip=zs[ngrip.i,ngrip.j,]
#best.dzs=zs.ngrip-zs.ngrip[pd]
#
##lines(time, best.dzs, col="#D71D5B", lty=2, lwd=2)
#
#axis(2,at=at.y,labels=at.y,col="grey30",col.axis="grey30",las=2, cex.axis=0.75)
#mtext("NGRIP Surface elevation anomaly (m)",side=2,line=3,cex=0.75, col="grey30")
#axis(1,at=at.x,labels=label.x,col="grey30",col.axis="grey30",las=1, cex.axis=0.75)
#mtext("Time (kyr ago)",side=1,line=2.5,cex=0.75, col="grey30")
#
## add GrIS
#file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
#nc=nc_open(file)
##u.obs=ncvar_get(nc,"uxy_srf")
#zs=raster(file, varname="z_srf")
#zb=raster(file, varname="z_bed")
#zb[zb<=0]=NA
#zs[zs<=0]=NA
#xc=ncvar_get(nc, "xc")
#yc=ncvar_get(nc, "yc")
#nc_close(nc)
#
#par(fig=c(0.02,0.28,0.48,1),new=T)
#image(zs, col="grey90", xaxt= "n", yaxt= "n", xlab="", ylab="",frame.plot=F)
#contour(zs, drawlabel=F, add=T, nlevels=5, col="grey40", lwd=0.5)
#points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1)
#
#dev.off()
#
#
#########################################################################################################################
# Plot dS + inset GrIS
plot.out=paste(work.fldr,"/Figures/fig5a-dS-NGRIP.png", sep="")
png(plot.out, width = 7, height = 5, units = "in", res=100)
#dev.new(width=7, height=5, units="in")
par(mar = c(0,0,0,0), oma = c(3.3,3.5,1,1),mgp=c(0,0.5,0))
#par(mar = c(1,2,1,1)) #Add a space 

#col=rev(plasma(nr.best))
#mycol=col[as.numeric(cut(b.best[,10],breaks = nr.best))]

# Plot dS
par(fig=c(0,1,0,1))
at.x=seq(-12,0,by=2)
label.x=seq(12,0,by=-2)
at.y=seq(-50,250, by=50)
xlim=c(-12,0)
ylim=c(-50,250)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)

polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("#C5E7D4",50))
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#B7E1C9",80))#78c6a3",50))
lines(t.ds.l13/1000, ds.ngrip1.l13, col=alpha.col("#7EC99E",30), lty=3,lwd=5)
lines(t.ds.l13/1000, ds.ngrip2.l13, col=alpha.col("#7EC99E",30), lty=3, lwd=5)
lines(t.ds.v09/1000, ds.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col=alpha.col("#7EC99E",50), lty=1, lwd=5)
#57CC99

# discarded runs
#for(i in nr.worst:1){
#  #file = paste(out.fldr,"/itmb.",itmb.worst[i],".itmc.",itmc.worst[i],".btq.",q.worst[i],".cbz0.",z0.worst[i],".cfngs.",cfngs.worst[i],".nffdlt.",neff.worst[i],"/yelmo2D.nc", sep="")
#  file=paste(out.fldr,"/",b.worst[i,1]-1,"/yelmo2D.nc",sep="")
#
#  # does the file exist?
#  io=file.exists(file)
#  if(io==F){
#  cat("Simulation ",file,"did not start \n")
#    #nom[i,]=NA
#    next
#  }
#  
#  nc = nc_open(file)
#  xc=ncvar_get(nc, "xc")
#  yc= ncvar_get(nc, "yc")
#  zs=ncvar_get(nc,"z_srf")
#  basins =ncvar_get(nc,"basins")
#  time=ncvar_get(nc,"time")
#  time=time/1000
#  #time=time-1.95
#  pd=length(time)
#  nc_close(nc)
#  zs.ngrip=zs[ngrip.i,ngrip.j,]
#  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
#  lines(time,  dzs.ngrip, col=alpha.col("grey20",100), lty=1, lwd=0.2)
#  #lines(time,  dzs.ngrip, col=alpha.col("#1E855A",100), lty=1, lwd=0.2)
#}

# normalization between 50 and 100 for colour opacity
b=100  # max opacity
a=5 # min opacity

# best runs
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  zs=ncvar_get(nc,"z_srf")
  basins =ncvar_get(nc,"basins")
  time=ncvar_get(nc,"time")
  time=time/1000
  #time=time-1.95
  pd=length(time)
  nc_close(nc)
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  #lines(time,  dzs.ngrip, col=alpha.col(mycol[i],ss.norm[i]*100), lty=1, lwd=2)
  lines(time,  dzs.ngrip, col=alpha.col(mycol[i],(b-a)*(ss.norm[i]-min(ss.norm))/(max(ss.norm)-min(ss.norm))+a ), lty=1, lwd=2)
}

# mean + 2*sd
lines(time,  mean.dzs, col=alpha.col("grey20",100), lty=1, lwd=3)
#lines(time,  mean.dzs+2*sd.dzs, col=alpha.col("grey20",100), lty=2, lwd=2)
#lines(time,  mean.dzs-2*sd.dzs, col=alpha.col("grey20",100), lty=2, lwd=2)
lines(time,  mean.dzs+1*sd.dzs, col=alpha.col("grey20",100), lty=2, lwd=3)
lines(time,  mean.dzs-1*sd.dzs, col=alpha.col("grey20",100), lty=2, lwd=3)

file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc",sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
zs=ncvar_get(nc,"z_srf")
nc_close(nc)
zs.ngrip=zs[ngrip.i,ngrip.j,]
dzs.ngrip=zs.ngrip-zs.ngrip[pd]
#lines(time, dzs.ngrip, col="#1982C4", lty=1, lwd=3)

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=0, cex.axis=1.3, tck=-0.02)
mtext("NGRIP Surface elevation anomaly (m)",side=2,line=2.4,cex=1.3)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1.3, tck=-0.015)
mtext("Time (kyr ago)",side=1,line=2.2,cex=1.3)

file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-16KM/GRL-16KM_TOPO-M17.nc"
nc=nc_open(file)
#u.obs=ncvar_get(nc,"uxy_srf")
zs=raster(file, varname="z_srf")
zb=raster(file, varname="z_bed")
zb[zb<=0]=NA
zs[zs<1]=NA
xc=ncvar_get(nc, "xc")
yc=ncvar_get(nc, "yc")
nc_close(nc)

#par(fig=c(0.02,0.30,0.48,1),new=T)
#image(zs, col="grey90", xaxt= "n", yaxt= "n", xlab="", ylab="",frame.plot=F)
#contour(zs, drawlabel=F, add=T, nlevels=5, col="grey40")
#points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1)

#par(fig=c(0.6,1,0.42,0.72),new=T)
par(fig=c(0.6,1,0.42,0.82),new=T)
legend_image <- as.raster(matrix(col, nrow=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')#, main = 'Thinning (m)', cex.main=0.8)
text(0.95,0.8,"Thinning (m)", cex=1.2)
labels=seq(round(min(b.best[1:nr.best,19]),digits=0),round(max(b.best[1:nr.best,19]),digits=0),l=5)
#text(x = seq(0.4,1.7,l=5), y=0.85, labels = round(labels), cex=0.6)
#rasterImage(legend_image, 0.3,0.45,1.7,0.7)
text(x = seq(0.4,1.7,l=5), y=0.65, labels = round(labels), cex=1)
rasterImage(legend_image, 0.3,0.45,1.7,0.6)
dev.off()
