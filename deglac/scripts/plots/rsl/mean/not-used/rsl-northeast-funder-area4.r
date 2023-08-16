library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)
#library(signal)

source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")

# North Funder et al., 2011 (Area 4 - 81.5°N) Prinsesse I ngeborg Halvø
loc=findij(81.537614, -13.791761)

t=c(-6840,-5509,-5039,-5022,-5006,-4740,-4485,-4476,-4258,-4192,-3772,-3123,-2805,-2760,-2541,-1819,-1810,-1781,-905,-815,-765,-756,-754,-734,-661,-646,-607,-567,-129,-126,-124,-100,-79,-67,-62,-10892,-10331,-9856,-10011,-9524,-9300,-8728,-8608,-7908,-7851,-7327,-6719,-6097,-3332)
t.err=c(142,220,237,288,262,540,445,305,266,238,199,207,183,160,185,219,123,206,274,115,144,299,151,188,266,261,128,80,142,155,157,183,189,211,400,372,414,311,385,381,166,279,238,110,114,151,265,188,113)
rsl=c(2,2,1,1,1,6,12,2,10,2,2,2,1,2,15,2,2,2,1,1,0,1,1,1,2,0,2,2,1,0,0,0,1,2,2,54,52,49,49,58,42,39,34,20,25,25,20,3,1)

xlim=c(-15000,0)
ylim=c(-50,200)
col.data="black"
col.mod=c("#26547c","#ef476f","#ffd166")

plot.out=paste(work.fldr,"/Figures/rsl/rsl-northeast-funder-area4.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new()
plot(NULL, xlim=xlim, ylim=ylim, ylab="RSL (m)", xlab="Time (yr)", main="Northeast - 81.5°N (Funder 2011)")

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
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area4.txt",sep=""), skip=1)
mean.rsl.b18=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area4.txt",sep=""), skip=1)
sd.rsl.b18=unlist(sd.rsl, use.names=FALSE)

time=read.table(paste(work.fldr,"/scripts/score-weighting/time.txt",sep=""), skip=1)
time=unlist(time, use.names=FALSE)
ty=time*1000

# load mean+ sd values B20
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac/B20")
mean.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/mean-rsl-northeast-funder-area4.txt",sep=""), skip=1)
mean.rsl.b20=unlist(mean.rsl, use.names=FALSE)
sd.rsl=read.table(paste(work.fldr,"/scripts/score-weighting/rsl/sd-rsl-northeast-funder-area4.txt",sep=""), skip=1)
sd.rsl.b20=unlist(sd.rsl, use.names=FALSE)

# B20
#polygon(c(ty, rev(ty)), c(mean.rsl.b20-2*sd.rsl.b20, rev(mean.rsl.b20+2*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b20-1*sd.rsl.b20, rev(mean.rsl.b20+1*sd.rsl.b20)), border = NA, col=alpha.col("#ffbd00",50))
lines(ty, mean.rsl.b20, col="#ffbd00", lwd=3, lty=1)

# B18
#polygon(c(ty, rev(ty)), c(mean.rsl.b18-2*sd.rsl.b18, rev(mean.rsl.b18+2*sd.rsl.b18)), border = NA, col=alpha.col("#390099",30))
polygon(c(ty, rev(ty)), c(mean.rsl.b18-1*sd.rsl.b18, rev(mean.rsl.b18+1*sd.rsl.b18)), border = NA, col=alpha.col("#390099",50))
lines(ty, mean.rsl.b18, col="#390099", lwd=3, lty=1)


for(i in 1:length(t)){
  arrows(x0=(t[i]-t.err[i]),y0=rsl[i],x1=(t[i]+t.err[i]),y1=rsl[i], code=3, angle=90, length=0.05, col=col.data, lwd=1)
}

dev.off()
