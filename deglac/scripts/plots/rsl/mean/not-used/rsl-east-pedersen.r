library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")

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

plot.out=paste(work.fldr,"/Figures/rsl/rsl-east-pedersen.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="EAST - Pedersen 2011")

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
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-east-pedersen.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-east-pedersen.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl, use.names=FALSE)

time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time*1000

# load mean+ sd values B20
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac/B20")
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-east-pedersen.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-east-pedersen.txt",sep=""), skip=1)
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
