library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
library(signal)

source("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij-8km.r")

# North East Bennike&WEidick 2001 (Midgardsormen)
loc=findij(79.65,-21.13)
t.min=c(-6620,-6670,-6410,-6420,-6400,-5200,-6400,-5590,-5440,-6890,-7150,-7400,-7170,-6290,-6300,-5980,-5660)
t.max=c(-6410,-6460,-6290,-6290,-6210,-4870,-6320,-5320,-5070,-6740,-6890,-7280,-6910,-6100,-6180,-5760,-5490)
rsl=c(11,11,17,11,17,15,17.5,10,12,20,30,12,NA,10,10,10,10)

xlim=c(-15000,0)
ylim=c(-50,400)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
#plot.out=paste("/home/itabone/work/v1.23/lhs/rsl-ne-blaso.png", sep="")
#png(plot.out, width=5, height=5, units="in", res=100)

plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - Midgardsormen (Bennike, Weidick 2001)")

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test2")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test2")
nbands=42 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
mycol=col[as.numeric(cut(b.best[,7],breaks = nr))] # ngrip.max

for(i in 1:nr){
  file = paste(out.fldr,"/btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
    arrows(x0=t.min[i],y0=rsl[i],x1=t.max[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
    #arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
  }

#dev.off()
####################################################################################################################
