library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")

# West Arveprinsens Ejland (Long 1999)
loc=findij(69.784934, -51.258730)
t.min=c(-10181,-8980,-8151,-7317,-5310,-4503,-1407)
t.max=c(-9563,-8387,-7511,-7156,-4980,-3721,-1066)
rsl=c(68.45,45.91,24.36,17.11,12.97,-2.75,-2.8)
rsl.err=c(0.46,0.38,0.78,0.25,0.93,0.37,0.37)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="WEST -  Arveprinsens Ejland")

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
nr=1

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

# plot data
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
    arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
}


file="/home/itabone/models/yelmo/v1.23/yelmox/output/test/B18-tas-B20-highP/kppgrz.5.fpn.0.02.cfn.1.0/yelmo2D.nc"
  nc = nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, col=col.mod[3],type="l", xlim=c(-15000,0))

