library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-marine-cores.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-cosmogenic-exposures.r")

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

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


###########################################################

### map

topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
bed=raster(topo, varname="z_bed")
surf=raster(topo, varname="z_srf")
mask=raster(topo, varname="mask")
#nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")


#file="/home/titan/gwgk/gwgk005h/work/79N/data/topography/BedMachineGreenland-2021-04-20-NE.nc"
#bed=raster(file, varname="bed")
#thick=raster(file, varname="thickness")
#surf=raster(file, varname="surface")
#mask=raster(file, varname="mask")


dev.new(width =4 , height = 4.5, units="in")
plot(bed)
contour(lat2D, nlevels=3, drawlabels=F, col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)


