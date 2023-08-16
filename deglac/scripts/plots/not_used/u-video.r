# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test7")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test7")
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


####################################################################################################################### Velocity fields

# mask areas outside greenland for grounding line 
topo.regions="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-16KM/GRL-16KM_REGIONS.nc"
region=raster(topo.regions, varname="mask")  # Ellesmere 1.11

topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-16KM/GRL-16KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
mask=raster(topo, varname="mask")
#nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")
#x2D=ncvar_get(nc,"x2D")
#y2D=ncvar_get(nc,"y2D")
#xc=ncvar_get(nc,"xc")
#yc=ncvar_get(nc,"yc")
#nc_close(nc)
# mask h over observed ice cover
htopo[htopo==0]=NA
# mask ellesmere island
htopo[region < 1.12 & region > 1.10]=NA

bands=c(11,12,27,37,47,57,67,87)

colors=colorRampPalette(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
                          "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
                          "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
                          "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
                          "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
                          "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
                          "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                          "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
                          "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
                          "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
                          "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
                          "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
                          "#AF0000", "#9F0000", "#8F0000", "#800000"))

colors=colorRampPalette(rev(c("#2E0F17","#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789","#1d3557")))

# plot
j=1
for(i in 1:75){
  plot.out=paste(work.fldr,"/Figures/video/u/u-video-",i,".png", sep="")
  png(plot.out, width = 6, height = 7, units = "in", res=100)
  #dev.new(width=6, height=7, units="in")
  #par(mfrow=c(2,4), mar=c(0.5,0.5,1.5,0.5), oma=c(1,1,1,5))
  #let=c("a","b","c","d","e","f","g","h")
  file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time = ncvar_get(nc,"time")
  nc_close(nc)
  #nc = nc_open(file)
  u = raster(file, varname="uxy_s", band=12+j)
  u[u==0.000] = NA
  u[u <0.1] = 0.1
  h = raster(file, varname="H_ice", band=12+j)
  h[h <= 0] = NA
  u[is.na(h)==T]=NA
  mask =raster(file, varname="mask_bed", band=12+j)
  u[mask!=2 & mask!=3]=NA
  bed = raster(file, varname="z_bed", band=12+j)
  bed[bed < 0]=NA

  if(i==nbands){
    plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, font.main=1, main=paste("Present day",sep=""), cex.main=1.5,box=F)
  }else{
    plot(bed, col=c("grey70","grey75","grey80","grey85","grey90"), legend=F, axes=F, font.main=1, main=paste(-time[12+j]/1000," kyr ago",sep=""), cex.main=1.5, box=F)
  }
  plot(log10(u), add=T, col=colors(100), zlim=c(log10(0.1), log10(1000)), legend=F, axes=F, box=F)
  contour(lat2D, nlevels=3, drawlabels=F, col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T) 
  contour(lon2D, nlevels=4, drawlabels=F,col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
  text(-460,-620,"80°",col="grey20")
  text(-160,-620,"60°",col="grey20")
  text(80,-620,"40°", col="grey20")
  text(340,-620,"20°", col="grey20")
  text(700,-620,"0° W", col="grey20")
  text(-640,-3322,"60° N",srt=-10, col="grey20")
  text(-640,-2150,"70° N",srt=-15, col="grey20")
  text(-640,-950,"80° N",srt=-25, col="grey20")
  #text(-680,-700, let[j], cex=2, font=2)  
  points(xc[ngrip.i],yc[ngrip.j],col="black",bg="black",pch=23, cex=1.5)  

  gris.range <- c(0.1, 1000)
   labels=c(0.1,1,10,100,1000)
 plot(log10(u), horizontal=F, legend.only=T, legend.shrink = 0.5, legend.width = 1.5, zlim=c(log10(0.1),log10(1000)),
      col=colors(100),
      axis.args=list(at=c(log10(0.1),log10(1),log10(10),log10(100),log10(1000)),
                     labels=labels,
                     cex.axis=1),
      legend.args=list(text='U (m/yr)', side=3, font=2, line=0.5, cex=1)
  )
  j = j+1
  dev.off()
}

################# legend ###########################
#plot.out=paste(work.fldr,"/Figures/U-hol-legend.png", sep="")
#png(plot.out, width = 5, height = 5, units = "in", res=100)
##dev.new(width=5, height=5, units="in")
# gris.range <- c(0.1, 1000)
#   labels=c(0.1,1,10,100,1000)
# plot(log10(u), horizontal=F, legend.only=T, legend.shrink = 0.5, legend.width = 1.5, zlim=c(log10(0.1),log10(1000)),
#      col=colors(100),
#      axis.args=list(at=c(log10(0.1),log10(1),log10(10),log10(100),log10(1000)),
#                     labels=labels,
#                     cex.axis=1),
#      #legend.args=list(text='U (m/yr)', side=3, font=2, line=0.5, cex=0.8)
#  )

#dev.off()
