library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij-8km.r")

# North Funder et al. 2011 (Area 5 - 80.5°N) Holm Land
loc=findij(80.421813, -15.791230)
t=c(-3654,-12687,-10887,-10349,-10231,-10195,-10074,-9525,-9322,-9264,-9251,-8565,-7911,-7301,-5053)
t.err=c(177,265,288,148,325,280,158,240,196,225,219,152,117,121,212)
rsl=c(10,22,34,62,65,30,54,50,15,22,46,18,6,13,4)

t.plant=c(-8142,-7416)
t.err.plant=c(180,152)
rsl.plant=c(39,26)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
plot.out=paste("/home/itabone/work/v1.23/lhs/B18/8km/test3/Figures/rsl-ne-area5.png", sep="")
png(plot.out, width=5, height=5, units="in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 80.5°N (Funder 2011)")

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
for(i in 1:length(t)){
  arrows(x0=(t.plant[i]-t.err.plant[i]),y0=rsl.plant[i],x1=(t.plant[i]+t.err.plant[i]),y1=rsl.plant[i], code=3, angle=90, length=0.05, col="grey50", lwd=1)
}

dev.off()
