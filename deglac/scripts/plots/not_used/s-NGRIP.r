# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr))
mycol=col[as.numeric(cut(b.best[,9],breaks = nr))]

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

###########################################################################################

# Plot S vs time
plot.out=paste(work.fldr,"/Figures/S-vs-time.png", sep="")
png(plot.out, width = 7, height = 5, units = "in", res=100)
#dev.new(width = 7, height = 5, units = "in")
par(mar=c(3,3,1,1), oma=c(2,2,0,0))
at.x=seq(-15,0,by=2.5)
label.x=seq(15,0,by=-2.5)
at.y=seq(2850,3150, by=50)
xlim=c(-15,0)
ylim=c(2850,3150)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
lines(t.ds.v09/1000, s.ngrip.v09, type="l", col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)), border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, s.ngrip1.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, s.ngrip2.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]

  lines(time,  zs.ngrip, col=alpha.col(mycol[i],100), lty=1, lwd=1.5)
}

axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=2, cex.axis=1)
mtext("Surface elevation (m)",side=2,line=3.5,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)",side=1,line=3,cex=1)
dev.off()

