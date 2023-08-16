library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
library(fields)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,19],breaks = nr.best))]
plot(b.best[,2],b.best[,17],pch=21, cex=1, col="black", bg=mycol)



################
#col=rev(plasma(nr.best))
col=rev(colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr))
#col=rev(tim.colors(nr.best))
b[,17][b[,17]<=0]=1e-20
mycol=col[as.numeric(cut(log(b[,17]),breaks = nr))] # mycol is wrt to score

b.small=b[,2:16]
pairs(b.small, pch=20,cex=0.5,col=mycol, upper.panel=NULL)

####### only best runs ######
col=rev(colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr.best))
b.best[,17][b.best[,17]<=0]=1e-20
mycol=col[as.numeric(cut(log(b.best[,17]),breaks = nr.best))] # mycol is wrt to score
#
b.small=b.best[,2:16]
pairs(b.small, pch=20,cex=0.5,col=mycol, upper.panel=NULL)
