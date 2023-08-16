library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

# North Funder et al., 2011 (Area 2 - 83°N - Cape Ole)
loc=findij(83.374692, -28.079069)
t=c(-7301,-10124,-9938,-9898,-9718,-9716,-9633,-9225,-8786,-8375,-8313,-8068,-6804,-6399,-4281)
t.err=c(123,162,229,267,169,171,111,190,195,163,112,98,104,97,129)
rsl=c(16,12,24,30,13,26,27,25,8,10,4,13,4,2,2)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 83°N (Funder 2011)")

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

#dev.off()
