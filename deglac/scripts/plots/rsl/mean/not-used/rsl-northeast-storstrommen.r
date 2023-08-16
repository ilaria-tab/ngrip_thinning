library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")

# North East Weidick 1996 (Storsrtrommen)
loc=findij(77.17,-22)
t.max=c(-1327,-3090,-3585,-3620,-4231,-5244,-5264,-5410,-5583)
t.min=c(-1091,-2702,-3130,-3334,-3872,-4714,-4830,-4959,-5135)
rsl=c(NA,13,140,NA,NA,150,150,NA,140)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")

plot.out=paste(work.fldr,"/Figures/rsl/rsl-northeast-storstrommen.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - Storstrommen (Weidick 1996)")

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
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-northeast-storstrommen.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-northeast-storstrommen.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl, use.names=FALSE)

time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time*1000

# load mean+ sd values B20
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac/B20")
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-northeast-storstrommen.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-northeast-storstrommen.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl, use.names=FALSE)

# B20
#polygon(c(ty, rev(ty)), c(mean.rsl.b20-2*sd.rsl.b20, rev(mean.rsl.b20+2*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# B18
#polygon(c(ty, rev(ty)), c(mean.rsl.b18-2*sd.rsl.b18, rev(mean.rsl.b18+2*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",50))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)


for(i in 1:length(t.max)){
   arrows(x0=t.max[i],y0=rsl[i],x1=t.min[i],y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

dev.off()
