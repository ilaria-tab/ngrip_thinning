library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test9_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test9_deglac")
nbands=87 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,13],breaks = nr.best))]

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

#######################################################################
#plot.out=paste(work.fldr,"/Figures/ne-tocn.png", sep="")
#png(plot.out, width = 6, height = 7, units = "in", res=100)
dev.new(width=5, height=5, units="in")
file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
tshlf=ncvar_get(nc,"tf_shlf")
nc_close(nc)

plot(time/1000, tshlf[156,290,],xlab="Time (kyr ago)",ylab="Tocn [°C]",type="l",lwd=5, ylim=c(-2,3),col="#ff006e", xlim=c(-15,0))
lines(time/1000, tshlf[127,322,], col="#8338ec", lwd=2)
lines(time/1000, tshlf[193,255,], col="#3a90ff", lwd=2)
grid()

