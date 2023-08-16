library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij-8km.r")

# North Funder et al., 2011 (Area 3 - 82.5°N - Herlufsholm Strand - Kap Clarence Wyckoff)
loc=findij(82.672789, -23.147292)
t=c(-6112,-5839,-5790,-4238,-9574,-9347,-9297,-9297,-9286,-13397,-8902,-8835,-8835,-8790,-8287,-7810,-7364,-7091,-6971,-5794,-5748,-3416,-3178,-3036)
t.err=c(193,145,135,388,103,254,168,171,193,293,214,187,186,31,98,560,95,150,384,128,144,141,205,166)
rsl=c(20,16,16,15,43,32,32,53,48,45,45,34,34,34,32,12,16,32,9,16,16,4,1,3)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
plot.out=paste("/home/itabone/work/v1.23/lhs/B18/8km/test3/Figures/rsl-ne-area3.png", sep="")
png(plot.out, width=5, height=5, units="in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 82.5°N (Funder 2011)")

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test3")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test3")
nbands=42 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
mycol=col[as.numeric(cut(b.best[,9],breaks = nr))]
for(i in 1:nr){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
