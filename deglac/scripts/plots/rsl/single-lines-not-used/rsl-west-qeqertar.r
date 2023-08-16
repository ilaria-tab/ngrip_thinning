library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")

# West Qeqertarsuatsiaq (Long and Roberts 2003)
loc=findij(68.462992, -52.909880)
t.min=c(-11041,-10547,-10668,-9272,-7647,-5259,-4567)
t.max=c(-10239,-10218,-10241,-8997,-7431,-4852,-4301)
rsl=c(NA,NA,97.09,50.7,31.18,7.41,1.89)
rsl.err=c(NA,NA,0.56,0.73,0.54,0.5,0.6)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="WEST - Qeqertarsuatsiaq")

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
#  lines(time, rsl.y, col=col.mod[1],type="l", xlim=c(-15000,0))
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


# Plot  data
for(i in 1:length(t.max)){
  arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
  arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
}







