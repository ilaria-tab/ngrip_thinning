library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

# East Pedersen 2011 DOI: 10.1002/jqs.1449
loc=findij(74.190652, -20.558731)
t=c(-6600,-7600,-7500,-4300,-5500,-4200,-580,-390,-180,-390,-2210)
t.err=c(300,400,600,300,400,300,20,40,40,40,50)
rsl=c(5.56,8.04,10.17,1.96,2.82,2.28,-0.28,-0.08,-0.14,-0.21,-0.53)
rsl.err=rep(0.22,13)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="EAST - Pedersen 2011")

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

# col alpha
alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}

#Plot data
for(i in 1:length(t)){
   arrows(x0=t[i]-t.err[i],y0=rsl[i],x1=t[i]+t.err[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
   arrows(x0=t[i],y0=rsl[i]-rsl.err[i],x1=t[i],y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
}

for(i in 1:nr.best){
  #file = paste(out.fldr,"/btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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


#Plot data
for(i in 1:length(t)){
   arrows(x0=t[i]-t.err[i],y0=rsl[i],x1=t[i]+t.err[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
   arrows(x0=t[i],y0=rsl[i]-rsl.err[i],x1=t[i],y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col=col.data, lwd=1)
}


