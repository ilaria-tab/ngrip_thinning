library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

###############################################################################################################

#file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
zs=raster(file, varname="z_srf")
zs[zs<0]=NA

file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa-negis.nc"
basin=raster(file, varname="basin")
#basin[basin!=9]=NA
negis=basin
#negis[negis!=9]=NA

topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
mask=raster(topo, varname="mask")
zb=raster(topo,varname="z_bed")
nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")
#x2D=ncvar_get(nc,"x2D")
#y2D=ncvar_get(nc,"y2D")
xc=ncvar_get(nc,"xc")
yc=ncvar_get(nc,"yc")
nc_close(nc)

file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17-paleonegis-1sd.nc"
nc=nc_open(file)
zb.negis=raster(file,varname="z_bed")
d=zb.negis-zb
d[d>=0]=NA


# NE transect-> 79N, Blaso
#xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698) # y coord
xs=c(200,744)
ys=c(-1354,-690)
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat.new=data.frame(x.line, y.line)

xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
ys=c(-1322,-698) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat.old=data.frame(x.line, y.line)

plot.out=paste(work.fldr,"/Figures/figS01-map-1D.png", sep="")
png(plot.out, width=4, height=5, units="in", res=100)
#dev.new(width=4, height=5, units="in")
par(mar=c(1,1,1,1))

plot(zs, col=c("#dee2e6"), box=F,axes=F,legend=F)
contour(zs, add=T, drawlabels=F, col="grey60")
#plot(negis, col=alpha.col("grey20",30), add=T,legend=F)
plot(d, add=T, zlim=c(-150,0), )
contour(negis, level=c(9),add=T, drawlabels=F,lwd=1,col="grey20")
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
  contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
points(xc[ngrip.i],yc[ngrip.j],bg="grey20",pch=23, cex=1,lwd=1)
#points(xc[egrip.i],yc[egrip.j],col="#6CC25B", pch=20, cex=2)
#points(xc[149],yc[305],bg="grey40", pch=21, cex=1, lwd=1)
#segments(x0=200, y0=-1354, x1=744 , y1=-690,col="orange", lty=1, lwd =2) # new segment
segments(x0=288, y0=-1322, x1=640 , y1=-698,col="magenta", lty=1, lwd =2) # old segment
#segments(x0=304, y0=-1322, x1=640 , y1=-698,col="magenta", lty=1, lwd =2) # old segment
#segments(x0=288, y0=-1322, x1=704, y1=-698, col="magenta", lty=1, lwd=2)
#segments(x0=272, y0=-1322, x1=704, y1=-698, col="magenta", lty=1, lwd=2)
points(xc[149],yc[305],bg="grey40", pch=21, cex=1, lwd=1)
text(0,-1700,"NGRIP",cex=0.8,font=2)
text(140,-1400,"NEGIS",cex=0.8,font=2)
text(650,-1000,"79N",cex=0.8,font=2)
#points(xc[156],yc[306],col="#1098F7", pch=20, cex=2)
dev.off()

