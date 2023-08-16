library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
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
alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}


#####################################################################################
#plot.out=paste(work.fldr,"/Figures/map-1D.png", sep="")
#png(plot.out, width=4, height=5, units="in", res=100)
dev.new(width=4, height=5, units="in")
par(mar=c(1,1,1,1))

#file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
zs=raster(file, varname="z_srf")
zs[zs<0]=NA

file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_BASINS-nasa-negis.nc"
basin=raster(file, varname="basin")
#basin[basin!=9]=NA
negis=basin
negis[negis!=9]=NA

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

plot(zs, col=c("#dee2e6"), box=F,axes=F,legend=F)
contour(zs, add=T, drawlabels=F, col="grey60")
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
  contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
#segments(x0=288, y0=-1322, x1=640 , y1=-698,col="magenta", lty=1, lwd =2) # old segment
text(0,-1500,"NGRIP",cex=0.8,font=2)
text(170,-1300,"NEGIS",cex=0.8,font=2)

# n-ne
#nne=findij(81.421896, -18.373036)
#points(nne[3],nne[4],bg="grey20",pch="1", cex=1,lwd=1)
blason=findij(79.62,-22.50)
points(blason[3],blason[4],bg="steelblue",pch="2", cex=1,lwd=1)
#blaso=findij(79.56,-22.38)
#points(blaso[3],blaso[4],bg="green",pch="3", cex=1,lwd=1)
fun1=findij(83.6,-31)
points(fun1[3],fun1[4],bg="magenta",pch="4", cex=1,lwd=1)
fun2=findij(83.374692, -28.079069)
points(fun2[3],fun2[4],bg="yellow",pch="5", cex=1,lwd=1)
fun3=findij(82.672789, -23.147292)
points(fun3[3],fun3[4],bg="blue",pch="6", cex=1,lwd=1)
fun4=findij(81.537614, -13.791761)
points(fun4[3],fun4[4],bg="orange",pch="7", cex=1,lwd=1)
fun5=findij(80.421813, -15.791230)
points(fun5[3],fun5[4],bg="pink",pch="8", cex=1,lwd=1)
hovn=findij(79.816667,-19.55)
points(hovn[3],hovn[4],bg="red",pch="9", cex=1,lwd=1)
#hovs=findij(79.75,-18.8)
#points(hovs[3],hovs[4],bg="lightblue",pch="$", cex=1,lwd=1)
ile=findij(77.83,-17.67)
points(ile[3],ile[4],bg="lightgreen",pch="%", cex=1,lwd=1)
#lam=findij(79.196315, -19.959059)
#points(lam[3],lam[4],bg="coral",pch="+", cex=1,lwd=1)
mid=findij(79.65,-21.13)
points(mid[3],mid[4],bg="darkorange",pch=23, cex=1,lwd=1)
sano=findij(78.116667,-20.3)
points(sano[3],sano[4],bg="gold",pch=22, cex=1,lwd=1)
#san=findij(77.983333,-21.45)
#points(san[3],san[4],bg="snow",pch=21, cex=1,lwd=1)
stor=findij(77.17,-22)
points(stor[3],stor[4],bg="salmon",pch=24, cex=1,lwd=1)
ehall=findij(71.454149, -24.438857)
points(ehall[3],ehall[4],bg="violet",pch=25, cex=1,lwd=1)
eped=findij(74.190652, -20.558731)
points(eped[3],eped[4],bg="steelblue",pch="x", cex=1,lwd=1)

#dev.off()

