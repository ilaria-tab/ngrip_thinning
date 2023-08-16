library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")

# South Nanortalik Sparrenbom 2006b & Bennike 2002
loc=findij(60.234804, -45.013219)
t.min=c(-1165.5,-3000,-8140,-2740,-5350,-9300,-9630,-9500,-9600,-9530,-10650,-10240,-10250,-10600,-12050,-12150,-12650,-15150)
t.max=c(-738.5,-1200,-2300,-2340,-4650,-8600,-9470,-9030,-9000,-9130,-9550,-9750,-9700,-10210,-10550,-11150,-11250,-14150)
rsl=c(-0.26,-2.3,-2.96,-6.01,-7.46,-2.96,-6.01,-7.46,-3.81,-1.6,3.49,1.49,2.49,9.49,22.49,24.49,28.49,31.49)
rsl.err=c(0.6,0.5,1.1,0.6,1.1,1.1,0.6,1.1,0.5,0.6,0.5,2,2,2,2,2,2,2)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="SOUTH - Nanortalik (Sparrenbom 2006b)")

precip=c("low","mean","high")
# B18
for(i in precip){
  file=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/lhs/B18/B18-tas-B20-",i,"P-MAR3.5/yelmo2D.nc",sep="")
  nc=nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, col=col.mod[1],type="l", xlim=c(-15000,0))
#  for(i in 1:length(t.max)){
#    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
#    arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
#  }

}

# B20 S3T-high P, monthly, tuned for NGRIP and GRIP
out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/new/GHF-M18/negis-1sd/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B20//GHF-M18/negis-1sd/test5")
nbands=111 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,

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
  lines(time, rsl.y, col=col.mod[2],type="l", xlim=c(-15000,0))
#  for(i in 1:length(t.max)){
#    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
#    arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
#  }

}
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
    arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
  }



