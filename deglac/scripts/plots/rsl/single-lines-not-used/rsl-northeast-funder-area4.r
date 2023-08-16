library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij-8km.r")

# North Funder et al., 2011 (Area 4 - 81.5°N) Prinsesse I ngeborg Halvø
loc=findij(81.537614, -13.791761)

t=c(-6840,-5509,-5039,-5022,-5006,-4740,-4485,-4476,-4258,-4192,-3772,-3123,-2805,-2760,-2541,-1819,-1810,-1781,-905,-815,-765,-756,-754,-734,-661,-646,-607,
   -567,-129,-126,-124,-100,-79,-67,-62,-10892,-10331,-9856,-10011,-9524,-9300,-8728,-8608,-7908,-7851,-7327,-6719,-6097,-3332)
t.err=c(142,220,237,288,262,540,445,305,266,238,199,207,183,160,185,219,123,206,274,115,144,299,151,188,266,261,128,80,142,155,157,183,189,211,400,372,414,
   311,385,381,166,279,238,110,114,151,265,188,113)
rsl=c(2,2,1,1,1,6,12,2,10,2,2,2,1,2,15,2,2,2,1,1,0,1,1,1,2,0,2,2,1,0,0,0,1,2,2,54,52,49,49,58,42,39,34,20,25,25,20,3,1)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
plot.out=paste("/home/itabone/work/v1.23/lhs/B18/8km/test3/Figures/rsl-ne-area4.png", sep="")
png(plot.out, width=5, height=5, units="in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 81.5°N (Funder et al.,)")

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test3")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test3")
nbands=42 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
mycol=col[as.numeric(cut(b.best[,9],breaks = nr))]
for(i in 1:nr){
file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc",sep="")
nc = nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, col=mycol[i],lwd=1,type="l", xlim=c(-15000,0))
}


for(i in 1:length(t)){
  arrows(x0=(t[i]-t.err[i]),y0=rsl[i],x1=(t[i]+t.err[i]),y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

dev.off()
