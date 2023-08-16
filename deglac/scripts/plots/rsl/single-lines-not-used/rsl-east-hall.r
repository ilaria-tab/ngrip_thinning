library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")

# East Hall 2010 https://doi.org/10.1016/j.quascirev.2010.03.013
loc=findij(71.454149, -24.438857)
t=c(-10990,-9913,-9757,-9651,-9577,-9533,-10059,-10038,-9993,-9964,-9964,-9919,-9413,-9104,-8835,-9521,-10356,-10333,-10264,-10066,-10010,-9690,-9275,-10686,-10495,-10317,-10312,-9673,-9580,-10612,-10371,-10272,-8866,-8299)
t.err=c(193,210,201,149,100,82,146,132,163,177,177,202,99,115,135,105,123,121,128,153,155,156,138,140,121,111,140,143,110,95,140,102,166,101)
rsl=c(101,73,61,81,61,66,73,69,75,75,81,75,57,40,49,58,67,63,61,67,67,52,58,89,78,75,72,62,72,73,73,74,37,33)
rsl.err=rep(4,34)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")
dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="EAST - (Hall 2010)")

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

# Plot data
for(i in 1:length(t)){
  arrows(x0=t[i]-t.err[i],y0=rsl[i],x1=t[i]+t.err[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
  #arrows(x0=t[i],y0=rsl[i]-rsl.err[i],x1=t[i],y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
}
points(t,rsl,pch=18,col=col.data)


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


# Plot data
for(i in 1:length(t)){
  arrows(x0=t[i]-t.err[i],y0=rsl[i],x1=t[i]+t.err[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
  #arrows(x0=t[i],y0=rsl[i]-rsl.err[i],x1=t[i],y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
}
points(t,rsl,pch=18,col=col.data)



