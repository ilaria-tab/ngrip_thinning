library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")

# North East Weidick 1996 (Storsrtrommen)
loc=findij(77.17,-22)
t.max=c(-1327,-3090,-3585,-3620,-4231,-5244,-5264,-5410,-5583)
t.min=c(-1091,-2702,-3130,-3334,-3872,-4714,-4830,-4959,-5135)
rsl=c(NA,13,140,NA,NA,150,150,NA,140)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - Storstrommen (Weidick 1996)")

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/lhs/fric/test")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/fric/test")
nbands=111 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
mycol=col[as.numeric(cut(b.best[,6],breaks = nr))]

for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".enhshr.",enh[i],"/yelmo2D.nc", sep="")
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


for(i in 1:length(t.max)){
   arrows(x0=t.max[i],y0=rsl[i],x1=t.min[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
 }


