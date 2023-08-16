library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)

# fun find ij from lat lon
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")


# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

# load rsl data
source("/home/titan/gwgk/gwgk005h/work/ngrip/data/rsl/funder2011.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/data/rsl/bennike-weidick2001.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/data/rsl/pedersen2011.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/data/rsl/hall2010.r")

# location B18 results
out.fldr.b18=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr.b18=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")

# location B20 results
out.fldr.b20=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr.b20=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")

# simulation time
time=read.table(paste(work.fldr.b18,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time*1000


# plot

# plot parameters
xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
at.y=seq(-50,200,50)
at.x=seq(-15000,0,2500)
label.x=seq(15,0,-2.5)

plot.out=paste(work.fldr,"/Figures/fig6-rsl-all.png", sep="")
png(plot.out, width=15, height=7, units="in", res=100)
#dev.new(width=15, height=7, units="in")
par(mar=c(1,1,1,1))
layout(matrix(c(1,2,3,4,5,1,6,7,8,9,10,11,12,13,14), ncol=5, byrow=TRUE))
#layout(matrix(c(1,2,3,4,1,5,6,7,1,8,9,10,11,12,13,14), ncol=4, byrow=TRUE))
#layout.show(a)

par(oma=c(2, 2, 2, 2), mar=c(2, 2, 2, 2))

# GrIS map
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/GRL-8KM_TOPO-M17.nc"
zs=raster(file, varname="z_srf")
zs[zs<0]=NA

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

#plot(zs,box=F,axes=F,legend=T)
#extent(zs)=c(0, 964, -1000, -566)
contour(zs,axes=F, drawlabels=F, col="grey75")
contour(lat2D, nlevels=3, drawlabels=T, method="edge",vfont=c("sans serif", "bold"), col=alpha.col("grey20", 60), lwd=0.5,lty=2,add=T)
  contour(lon2D, nlevels=5, drawlabels=T, method="edge",vfont=c("sans serif", "bold"),col=alpha.col("grey20", 60), lwd=0.5,lty=2,add=T)
#text(-600,-500, "a", cex=1.5, font=2)
col.map=rainbow(12)
col.map=rep("steelblue",12)
#fun1=findij(83.6,-31)
points(loc.fun1[3],loc.fun1[4],bg=col.map[1],pch=23, cex=1,lwd=1)
text(loc.fun1[3]-50,loc.fun1[4], "a")
#fun2=findij(83.374692, -28.079069)
points(loc.fun2[3],loc.fun2[4],bg=col.map[2],pch=23, cex=1,lwd=1)
text(loc.fun2[3],loc.fun2[4]+70, "b")
#fun3=findij(82.672789, -23.147292)
points(loc.fun3[3],loc.fun3[4],bg=col.map[3],pch=23, cex=1,lwd=1)
text(loc.fun3[3],loc.fun3[4]+70, "c")
#fun4=findij(81.537614, -13.791761)
points(loc.fun4[3],loc.fun4[4],bg=col.map[4],pch=23, cex=1,lwd=1)
text(loc.fun4[3],loc.fun4[4]+70, "d")
#fun5=findij(80.421813, -15.791230)
points(loc.fun5[3],loc.fun5[4],bg=col.map[5],pch=23, cex=1,lwd=1)
text(loc.fun5[3]+70,loc.fun5[4], "e")
#hovn=findij(79.816667,-19.55)
points(loc.hovn[3],loc.hovn[4],bg=col.map[6],pch=23, cex=1,lwd=1)
text(loc.hovn[3]+70,loc.hovn[4], "f")
#mid=findij(79.65,-21.13)
points(loc.mid[3],loc.mid[4],bg=col.map[9],pch=23, cex=1,lwd=1)
text(loc.mid[3]+70,loc.mid[4]-30, "g")
#blason=findij(79.62,-22.50)
points(loc.blason[3],loc.blason[4],bg=col.map[7],pch=23, cex=1,lwd=1)
text(loc.blason[3]-60,loc.blason[4], "h")
#sano=findij(78.116667,-20.3)
points(loc.sano[3],loc.sano[4],bg=col.map[10],pch=23, cex=1,lwd=1)
text(loc.sano[3]-60,loc.sano[4], "i")
#ile=findij(77.83,-17.67)
points(loc.ile[3],loc.ile[4],bg=col.map[8],pch=23, cex=1,lwd=1)
text(loc.ile[3]+70,loc.ile[4], "j")
#eped=findij(74.190652, -20.558731)
points(loc.eped[3],loc.eped[4],bg=col.map[12],pch=23, cex=1,lwd=1)
text(loc.eped[3]+80,loc.eped[4], "k")
#ehall=findij(71.454149, -24.438857)
points(loc.ehall[3],loc.ehall[4],bg=col.map[11],pch=23, cex=1,lwd=1)
text(loc.ehall[3]+80,loc.ehall[4], "l")

#text(-600,-500, "a", cex=2, font=2)

#####################################################################################################
#####################################################################################################
###############     now plot rsl curves ###########################################
#################################################################################################


# funder1

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 83.5°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# mean + sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area1.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area1.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# mean + sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area1.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area1.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.fun1)){
  arrows(x0=(t.fun1[i]-t.err.fun1[i]),y0=rsl.fun1[i],x1=(t.fun1[i]+t.err.fun1[i]),y1=rsl.fun1[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}
for(i in 1:length(t.shell.moll)){
  arrows(x0=(t.shell.moll[i]-t.err.shell.moll[i]),y0=rsl.shell.moll[i],x1=(t.shell.moll[i]+t.err.shell.moll[i]),y1=rsl.shell.moll[i], code=3, angle=90, length=0.05, col="black", lwd=1)
}
for(i in 1:length(t.plant.moll)){
  arrows(x0=(t.plant.moll[i]-t.err.plant.moll[i]),y0=rsl.plant.moll[i],x1=(t.plant.moll[i]+t.err.plant.moll[i]),y1=rsl.plant.moll[i], code=3, angle=90, length=0.05, col="black", lwd=1)
}

text(-500,180, "a", cex=1.5, font=2)

##################################
# funder2:  North Funder et al., 2011 (Area 2 - 83°N - Cape Ole)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 83.3°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area2.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area2.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area2.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area2.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.fun2)){
  arrows(x0=(t.fun2[i]-t.err.fun2[i]),y0=rsl.fun2[i],x1=(t.fun2[i]+t.err.fun2[i]),y1=rsl.fun2[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "b", cex=1.5, font=2)

##################################
# funder 3:  North Funder et al., 2011 (Area 3 - 82.5°N - Herlufsholm Strand - Kap Clarence Wyckoff)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 82.6°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area3.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area3.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area3.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area3.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.fun3)){
  arrows(x0=(t.fun3[i]-t.err.fun3[i]),y0=rsl.fun3[i],x1=(t.fun3[i]+t.err.fun3[i]),y1=rsl.fun3[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "c", cex=1.5, font=2)


#############################################################
# funder4: North Funder et al., 2011 (Area 4 - 81.5°N) Prinsesse I ngeborg Halvø

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 81.5°N")
title("RSL (m) - 81.5°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area4.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area4.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area4.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area4.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.fun4)){
  arrows(x0=(t.fun4[i]-t.err.fun4[i]),y0=rsl.fun4[i],x1=(t.fun4[i]+t.err.fun4[i]),y1=rsl.fun4[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "d", cex=1.5, font=2)


#############################################
# funder 5: North Funder et al. 2011 (Area 5 - 80.5°N) Holm Land

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 80.4°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area5.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area5.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area5.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area5.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)


for(i in 1:length(t.fun5)){
  arrows(x0=(t.fun5[i]-t.err.fun5[i]),y0=rsl.fun5[i],x1=(t.fun5[i]+t.err.fun5[i]),y1=rsl.fun5[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}
for(i in 1:length(t.plant.moll)){
  arrows(x0=(t.plant.moll[i]-t.err.plant.moll[i]),y0=rsl.plant.moll[i],x1=(t.plant.moll[i]+t.err.plant.moll[i]),y1=rsl.plant.moll[i], code=3, angle=90, length=0.05, col="black", lwd=1)
}

text(-500,180, "e", cex=1.5, font=2)

###################################################
# hovsgard nord: North East Bennike&Weidick 2001 (Hovgaard Oer)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 79.8°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-hovgaard-north.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-hovgaard-north.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-hovgaard-north.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-hovgaard-north.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.max.hovn)){
    arrows(x0=t.min.hovn[i],y0=rsl.hovn[i],x1=t.max.hovn[i],y1=rsl.hovn[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "f", cex=1.5, font=2)

##############################################################
# blaso nord: North East Bennike&WEidick 2001 (Blaso NORTH)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 79.6°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-blaso-north.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-blaso-north.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-blaso-north.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-blaso-north.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.max.blason)){
    arrows(x0=t.min.blason[i],y0=rsl.blason[i],x1=t.max.blason[i],y1=rsl.blason[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "g", cex=1.5, font=2)

##########################################################
# midgaard:  North East Bennike&WEidick 2001 (Midgardsormen)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 79.6°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-midgardsormen.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-midgardsormen.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-midgardsormen.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-midgardsormen.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

for(i in 1:length(t.max.mid)){
    arrows(x0=t.min.mid[i],y0=rsl.mid[i],x1=t.max.mid[i],y1=rsl.mid[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "h", cex=1.5, font=2)

####################################################################
# temperature at 79N

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# dTann B18
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18-TAS-B20-highP.nc"
nc=nc_open(file)
dtas.b18=ncvar_get(nc, "tas")
t.b18=ncvar_get(nc, "time")
nc_close(nc)
dtas.b18.79=dtas.b18[149,305,,]
dtann.b18.79=c()
for(t in 1:length(t.b18)){
  dtann.b18.79[t]=mean(dtas.b18.79[1:12,t])
}
dtann.b18.79=dtann.b18.79-dtann.b18.79[2181] # shift the anomaly wrt to 1750 CE (t[2181])

# dTann B20
file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B20-S3T-highP-monthly.nc"
nc=nc_open(file)
dtas.b20=ncvar_get(nc, "tas")
t.b20=ncvar_get(nc, "time")
nc_close(nc)
dtas.b20.79=dtas.b20[149,305,,]
dtann.b20.79=c()
for(t in 1:length(t.b20)){
  dtann.b20.79[t]=mean(dtas.b20.79[1:12,t])
}
dtann.b20.79=dtann.b20.79-dtann.b20.79[397] # shift the anomaly wrt to 1750 CE (t[397])

at.y.temp=seq(-12.5,7.5, by=2.5)
xlim.temp=c(-15,0)
ylim.temp=c(-12.5,7.5)
at.x.temp=seq(-15,0, by=2.5)
lab.x.temp=seq(15,0, by=-2.5)
plot(xlim.temp, ylim.temp, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
grid()
polygon(c(t.tann.l17/1000, rev(t.tann.l17/1000)), c(tann.agas.max.l17, rev(tann.agas.min.l17)), border = NA, col=alpha.col("grey60",30))
lines(t.tann.l17/1000,tann.agas.l17, lwd=2, col="grey60")
lines(t.b18/1000, dtann.b18.79, lwd=1.5, col="#390099") #"#0496ff")
lines(t.b20/1000, dtann.b20.79, lwd=2, col="#ffbd00") #"#f72585")
axis(2,at=at.y.temp,labels=at.y.temp,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Tann anomaly at 79N (°C)",side=2,line=2.3,cex=0.8)
axis(1,at=at.x.temp,labels=lab.x.temp,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)", side=1, line=2.3, cex=0.8)
#text(-12.5,7, "b", cex=1.5, font=2)



######################################################################
# san offshore: North East Bennike & Weidick 2001 (Sanddalen)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 78°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)
mtext("Time (kyr ago)", side=1, line=2.3, cex=0.8)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-sanddalen-offshore.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-sanddalen-offshore.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-sanddalen-offshore.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-sanddalen-offshore.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)


for(i in 1:length(t.max.sano)){
    arrows(x0=t.min.sano[i],y0=rsl.sano[i],x1=t.max.sano[i],y1=rsl.sano[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "i", cex=1.5, font=2)


##################################
# ile de france: North East Bennike & Weidick 2001 (Ile de France)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 77.8°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
mtext("Time (kyr ago)", side=1, line=2.3, cex=0.8)
abline(h=0, lty=2)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-northeast-ilefrance.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-northeast-ilefrance.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-northeast-ilefrance.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-northeast-ilefrance.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.max.ile)){
    arrows(x0=t.min.ile[i],y0=rsl.ile[i],x1=t.max.ile[i],y1=rsl.ile[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "j", cex=1.5, font=2)


###################################################################################
# east pedersen: East Pedersen 2011 DOI: 10.1002/jqs.1449

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="1. RSL (m) - 83.5°N")
title("RSL (m) - 74.2°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)
mtext("Time (kyr ago)", side=1, line=2.3, cex=0.8)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-east-pedersen.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-east-pedersen.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-east-pedersen.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-east-pedersen.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# RSL data
for(i in 1:length(t.eped)){
  arrows(x0=t.eped[i]-t.err.eped[i],y0=rsl.eped[i],x1=t.eped[i]+t.err.eped[i],y1=rsl.eped[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "k", cex=1.5, font=2)

#################################################################################################
# East Hall 2010 https://doi.org/10.1016/j.quascirev.2010.03.013

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", cex.axis=1, main="RSL (m) - 83.5°N")
title("RSL (m) - 71.5°N")
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1)
abline(h=0, lty=2)
mtext("Time (kyr ago)", side=1, line=2.3, cex=0.8)

# load mean+ sd values B18
mean.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/mean-rsl-east-hall.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl.b18, use.names=FALSE)
sd.rsl.b18=read.table(paste(work.fldr.b18,"/scripts/score-weighting/rsl/sd-rsl-east-hall.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl.b18, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)

# load mean+ sd values B20
mean.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/mean-rsl-east-hall.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl.b20, use.names=FALSE)
sd.rsl.b20=read.table(paste(work.fldr.b20,"/scripts/score-weighting/rsl/sd-rsl-east-hall.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl.b20, use.names=FALSE)
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# Plot data
for(i in 1:length(t.ehall)){
  arrows(x0=t.ehall[i]-t.err.ehall[i],y0=rsl.ehall[i],x1=t.ehall[i]+t.err.ehall[i],y1=rsl.ehall[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

text(-500,180, "l", cex=1.5, font=2)


dev.off()
