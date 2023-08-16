
library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)
library(TeachingDemos)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")


# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,9],breaks = nr.best))]

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

rsl.lec=c(-80,-60,-40,-10)
rsl.lec1=c(-40,-30,-20,-5)
time.lec=c(-16000,-12000,-8000,-4000)


plot(time, rsl.y, xlim=c(-18000,0), type="l")
points(time.lec, rsl.lec, pch=21, col="red")
points(time.lec, rsl.lec1, pch=21, col="blue")
