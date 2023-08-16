library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")

# West Akulliit/Nuuk (Long and Roberts 2002)
loc=findij(68.662742, -51.225047)
t.min=c(-9817,-8647,-7786,-6651,-8376,-8604,-7998,-7563)
t.max=c(-9436,-8390,-7571,-6310,-8032,-8395,-7742,-7266)
rsl=c(NA,56.56,35.84,16.4,56,56.38,43.62,33.27)
rsl.err=c(NA,0.59,0.22,0.79,0.22,0.82,0.69,0.57)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="WEST - Akulliit/Nuuk (Long and Roberts 2002)")

#precip=c("low","mean","high")
# B18
#for(i in precip){
#  file=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/lhs/B18/B18-tas-B20-",i,"P-MAR3.5/yelmo2D.nc",sep="")
#  nc=nc_open(file)
#  sl=ncvar_get(nc,"z_sl")
#  zb=ncvar_get(nc,"z_bed")
#  time=ncvar_get(nc,"time")
#  nc_close(nc)
#
#  sl.w=sl[loc[1],loc[2],]
#  zb.w=zb[loc[1],loc[2],]
#
#  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
#  lines(time, rsl.y, xlim=c(-15000,0), col=col.mod[1])
#}

# B18 
out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/new/B18/test1")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/test1")
nbands=111 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
nr=10

for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".kppgrz.",kgrz[i],".btq.",q[i],".cbz0.",z0[i],".cfn.",cfn[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, col=col.mod[1],type="l", xlim=c(-15000,0))
}


# B20 S3T-high P, monthly, tuned for NGRIP and GRIP
out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/new/GHF-M18/negis-1sd/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B20/GHF-M18/negis-1sd/test5")
nbands=111 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
# col palette
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
#mycol=col[as.numeric(cut(b.best[,9],breaks = nr))]
#nr=10
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".kppgrz.",kgrz[i],".fpn.",fpn[i],".btq.",q[i],".cbz0.",z0[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  if(b.ord[i,9] > 100){
    lines(time, rsl.y, col=col.mod[2],type="l", xlim=c(-15000,0),lwd=2)
  }else{
    next
  }
}

# Plot data
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
    arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
}

