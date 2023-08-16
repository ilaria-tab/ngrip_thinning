# plot showing NEGIS north/south switch and negis basin

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90  # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

## col alpha
#alpha.col <- function(col,percent=50)
#{
#  c <- rep("#000000",length(col))
#  for (i in 1:length(col)) {
#    cc <- col2rgb(col[i])/255
#    c[i] <- rgb(t(cc),alpha=percent/100)
#  }
#
#  return(c)
#}

#########################################################################################################

# function NEGIS north/south switch


# load files
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
zs=raster(file, varname="z_srf")
zs[zs<0]=NA

file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa-negis-three-2.nc"
basin=raster(file, varname="basin_sub")
#basin[basin!=9]=NA
negis=basin
negis[negis<9.1]=NA

topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
mask=raster(topo, varname="mask")
nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")
#x2D=ncvar_get(nc,"x2D")
#y2D=ncvar_get(nc,"y2D")
xc=ncvar_get(nc,"xc")
yc=ncvar_get(nc,"yc")
nc_close(nc)
###########################################################

# plots
plot.out=paste(work.fldr,"/Figures/fig1a-negis-switch.png", sep="")
png(plot.out, width = 5, height = 4, units = "in", res=100)
#dev.new(width = 5, height = 4, units="in")
#par(mfrow=c(1,2))
par(mar=c(3,4,1,0),oma=c(1,3,0,1))

at.x=seq(-12,0,by=2)
label.x=seq(12,0,by=-2)
at.y=c(0,0.3,1)
xlim=c(-12,0)
ylim=c(0,1)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n",yaxt="n", cex.axis=1.5)
segments(-12,1,-8.01,1, col="#01E63A",lwd=6)
segments(-8.01,1,-8,0.3, col="#01E63A", lwd=6)
segments(-8,0.3,0,0, col="#01E63A", lwd=6)
segments(-12,0.3,-8,0.3, col="#ff006e",lwd=6)
segments(-8,0.3,0,1, col="#ff006e", lwd=6)
segments(-12,0.3,0,0.3, col="#ffbe0b", lwd=6, lty=3)
text(-4,0.73,"North",cex=1.2, font=3, srt=43)
text(-4,0.35,"Centre",cex=1.2, font=3, srt=0)
text(-4,0.085,"South",cex=1.2, font=3, srt=-20)
axis(2,at=at.y,labels=c(expression("f"[min]),expression("f"[mid]),"1"),las=3,col="black",col.axis="black",las=1, cex.axis=1.3)
axis(1,at=at.x, labels=label.x, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr ago)",side=1,line=2.75,cex=1.3)

dev.off()

# second plot
plot.out=paste(work.fldr,"/Figures/fig1b-negis-switch-map.png", sep="")
png(plot.out, width = 5, height = 4, units = "in", res=100)
#dev.new(width = 5, height = 4, units="in")
#par(mfrow=c(1,2))
par(mar=c(0,0,0,0))
zs[zs<=0]=NA
plot(zs, col=alpha.col(colorRampPalette(brewer.pal(9,"Blues"))(100),80), box=F,axes=F,legend=F)
contour(zs, add=T, drawlabels=F, col="grey60")
#negis.c=negis[negis<9.2]
#plot(negis, col=alpha.col("grey20",30), add=T,legend=F)
plot(negis, col=alpha.col(c("#ffbe0b","#01E63A","#ff006e"),90), add=T,legend=F)
points(xc[ngrip.i],yc[ngrip.j],bg="grey20",pch=23, cex=1,lwd=1)
#points(xc[egrip.i],yc[egrip.j],bg="grey20",pch=21, cex=1,lwd=1)

contour(lat2D, nlevels=3, drawlabels=F, col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
contour(lon2D, nlevels=4, drawlabels=F,col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
text(-460,-620,"80",col="grey40", cex=0.8)
  text(-160,-620,"60",col="grey40", cex=0.8)
  text(80,-620,"40", col="grey40", cex=0.8)
  text(340,-620,"20", col="grey40", cex=0.8)
  text(700,-620,"0 W", col="grey40", cex=0.8)
  text(-640,-3322,"60 N",srt=-10, col="grey40", cex=0.8)
  text(-640,-2150,"70 N",srt=-15, col="grey40", cex=0.8)
  text(-640,-950,"80 N",srt=-25, col="grey40", cex=0.8)

text(0,-1530,"NGRIP",cex=0.8,font=1)
#text(400,-1650,"EGRIP",cex=0.8,font=1)
text(180,-1080,"North",cex=0.8, font=3, srt=25)
text(700,-1080,"Centre",cex=0.8, font=3, srt=0)
text(400,-1500,"South",cex=0.8, font=3, srt=60)

dev.off()

