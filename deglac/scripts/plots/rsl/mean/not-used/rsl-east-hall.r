library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
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
plot.out=paste(work.fldr,"/Figures/rsl/rsl-east-hall.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="EAST - (Hall 2010)")

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


# load mean+ sd values B18
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-east-hall.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-east-hall.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl, use.names=FALSE)

time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time*1000

# load mean+ sd values B20
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac/B20")
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-east-hall.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-east-hall.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl, use.names=FALSE)

# B20
#polygon(c(ty, rev(ty)), c(mean.rsl.b20-2*sd.rsl.b20, rev(mean.rsl.b20+2*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# B18
#polygon(c(ty, rev(ty)), c(mean.rsl.b18-2*sd.rsl.b18, rev(mean.rsl.b18+2*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",50))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)



# Plot data
for(i in 1:length(t)){
  arrows(x0=t[i]-t.err[i],y0=rsl[i],x1=t[i]+t.err[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
  #arrows(x0=t[i],y0=rsl[i]-rsl.err[i],x1=t[i],y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
}
points(t,rsl,pch=18,col=col.data)


dev.off()
