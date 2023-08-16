library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

# North Funder et al 2011 & Möller 2010 (Area 1 - 83.5°N)
loc=findij(83.6,-31)
t=c(-7930,-7847,-7292,-7106,-5813,-5795,-4105,-3720,-3297,-3223,-2497,-1808,-1173,-596,-141,-133,-122,-81,-78,-50)
t.err=c(83,148,133,145,177,130,274,163,147,161,216,123,159,58,144,294,163,224,205,253)
rsl=c(15,5,9,12,6,5,2,7,4,1,1,1,0,0,4,0,0,0,2,0)

t.shell=c(-10087,-10010,-9716,-9712,-9680,-9570,-9547,-10087,-10824,-9016,-8765,-8618,-8547,-8182,-8159,-7433,-7399,-7295,-7290,-7137,-7099,-4635,-1046)
t.err.shell=c(164,224,172,175,203,92,78,160,256,232,213,136,136,149,138,120,106,123,126,131,152,183,548)
rsl.shell=c(23,30,22,32,1,12,23,23,39,11,1,8,18,2,12,8,8,2,4,6,2,0,3)

t.plant=c(-10285,-9716,-7399,-6324,-5142,-7290,-5945,-5111)
t.err.plant=c(365,172,106,115,168,126,53,180)
rsl.plant=c(107,22,8,21,19,4,12,26)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 83.5°N (Funder 2011)")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?
source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best,
nr.best=nrow(b.best)
ss.norm=b.best[,17]
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr.best)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))] # ngrip.max

for(i in 1:length(t)){
  arrows(x0=(t[i]-t.err[i]),y0=rsl[i],x1=(t[i]+t.err[i]),y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}
for(i in 1:length(t.shell)){
  arrows(x0=(t.shell[i]-t.err.shell[i]),y0=rsl.shell[i],x1=(t.shell[i]+t.err.shell[i]),y1=rsl.shell[i], code=3, angle=90, length=0.05, col="grey30", lwd=1)
}
for(i in 1:length(t.plant)){
  arrows(x0=(t.plant[i]-t.err.plant[i]),y0=rsl.plant[i],x1=(t.plant[i]+t.err.plant[i]),y1=rsl.plant[i], code=3, angle=90, length=0.05, col="grey50", lwd=1)
}


for(i in 1:nr.best){
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc = nc_open(file)
  sl=ncvar_get(nc,"z_sl")
  zb=ncvar_get(nc,"z_bed")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  sl.w=sl[loc[1],loc[2],]
  zb.w=zb[loc[1],loc[2],]

  rsl.y=sl.w-zb.w-(sl.w[length(time)]-zb.w[length(time)])
  lines(time, rsl.y, col=alpha.col(mycol[i],ss.norm[i]*100),lwd=1,type="l", xlim=c(-15000,0))
}


for(i in 1:length(t)){
  arrows(x0=(t[i]-t.err[i]),y0=rsl[i],x1=(t[i]+t.err[i]),y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}
for(i in 1:length(t.shell)){
  arrows(x0=(t.shell[i]-t.err.shell[i]),y0=rsl.shell[i],x1=(t.shell[i]+t.err.shell[i]),y1=rsl.shell[i], code=3, angle=90, length=0.05, col="grey30", lwd=1)
}
for(i in 1:length(t.plant)){
  arrows(x0=(t.plant[i]-t.err.plant[i]),y0=rsl.plant[i],x1=(t.plant[i]+t.err.plant[i]),y1=rsl.plant[i], code=3, angle=90, length=0.05, col="grey50", lwd=1)
}

