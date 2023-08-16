library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-marine-cores-4km.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-cosmogenic-exposures-4km.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-4km.r")

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

###########################################################
# cosmogenic exposure data
bourb.age=bourb.age.min+(bourb.age.max-bourb.age.min)/2
bourb.err=(bourb.age.max-bourb.age.min)/2
storoen.age=storoen.age.min+(storoen.age.max-storoen.age.min)/2
storoen.err=(storoen.age.max-storoen.age.min)/2
kapam.age=kapam.age.min+(kapam.age.max-kapam.age.min)/2
kapam.err=(kapam.age.max-kapam.age.min)/2
blaso.age=blaso.age.min+(blaso.age.max-blaso.age.min)/2
blaso.err=(blaso.age.max-blaso.age.min)/2
lambland1.age=lambland1.age.min+(lambland1.age.max-lambland1.age.min)/2
lambland1.err=(lambland1.age.max-lambland1.age.min)/2
lambland2.age=lambland2.age.min+(lambland2.age.max-lambland2.age.min)/2
lambland2.err=(lambland2.age.max-lambland2.age.min)/2
zi.age=zi.age.min+(zi.age.max-zi.age.min)/2
zi.err=(zi.age.max-zi.age.min)/2
sondrem.age=sondrem.age.min+(sondrem.age.max-sondrem.age.min)/2
sondrem.err=(sondrem.age.max-sondrem.age.min)/2
blochn.age=blochn.age.min+(blochn.age.max-blochn.age.min)/2
blochn.err=(blochn.age.max-blochn.age.min)/2


# marine sediments
ps100.deg.time=c(-10200,-9600)
g92.deg.time=c(-13400,-11200)
g39.deg.time=c(-14000,-13300)
ps100198.deg.time=c(-10900,-10000)

ps100.age=(-10200)+((-9600)-(-10200))/2
ps100.err=((-9600)-(-10200))/2
g92.age=(-13400)+((-11200)-(-13400))/2
g92.err=((-11200)-(-13400))/2
g39.age=(-14000)+((-13300)-(-14000))/2
g39.err=((-13300)-(-14000))/2
ps100198.age=(-10900)+((-10000)-(-10900))/2
ps100198.err=((-10000)-(-10900))/2

################################################################################################################
# B18

mis.ps100=read.table(paste(work.fldr,"/rmse-ps100.txt",sep=""), skip=1)
mis.g92=read.table(paste(work.fldr,"/rmse-g92.txt",sep=""), skip=1)
mis.g39=read.table(paste(work.fldr,"/rmse-g39.txt",sep=""), skip=1)
mis.ps100198=read.table(paste(work.fldr,"/rmse-ps100198.txt",sep=""), skip=1)

mis.ps100=as.vector(mis.ps100[,1])
mis.g92=as.vector(mis.g92[,1])
mis.g39=as.vector(mis.g39[,1])
mis.ps100198=as.vector(mis.ps100198[,1])

mis.kapam=read.table(paste(work.fldr,"/rmse-kapam.txt",sep=""), skip=1)
mis.storoen=read.table(paste(work.fldr,"/rmse-storoen.txt",sep=""), skip=1)
mis.bourb=read.table(paste(work.fldr,"/rmse-bourb.txt",sep=""), skip=1)
mis.blochn=read.table(paste(work.fldr,"/rmse-blochn.txt",sep=""), skip=1)
mis.sondrem=read.table(paste(work.fldr,"/rmse-sondrem.txt",sep=""), skip=1)
mis.zi=read.table(paste(work.fldr,"/rmse-zi.txt",sep=""), skip=1)
mis.lambland1=read.table(paste(work.fldr,"/rmse-lambland1.txt",sep=""), skip=1)
mis.lambland2=read.table(paste(work.fldr,"/rmse-lambland2.txt",sep=""), skip=1)
mis.blaso=read.table(paste(work.fldr,"/rmse-blaso.txt",sep=""), skip=1)

mis.kapam=as.vector(mis.kapam[,1])
mis.storoen=as.vector(mis.storoen[,1])
mis.bourb=as.vector(mis.bourb[,1])
mis.blochn=as.vector(mis.blochn[,1])
mis.sondrem=as.vector(mis.sondrem[,1])
mis.zi=as.vector(mis.zi[,1])
mis.lambland1=as.vector(mis.lambland1[,1])
mis.lambland2=as.vector(mis.lambland2[,1])
mis.blaso=as.vector(mis.blaso[,1])


b18.ps100=mis.ps100
b18.g92=mis.g92
b18.g39=mis.g39
b18.ps100198=mis.ps100198
b18.kapam=mis.kapam
b18.storoen=mis.storoen
b18.bourb=mis.bourb
b18.blochn=mis.blochn
b18.sondrem=mis.sondrem
b18.zi=mis.zi
b18.lambland1=mis.lambland1
b18.lambland2=mis.lambland2
b18.blaso=mis.blaso


# weighted mean
ss.norm=b.best[,17]
# median(b18.g92[b.best[,1]]

mean.b18.ps100=sum(ss.norm*b18.ps100[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.ps100=sqrt(sum(ss.norm*(b18.ps100[b.best[,1]]-mean.b18.ps100)^2)/(w*sum(ss.norm)))

mean.b18.g92=sum(ss.norm*b18.g92[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.g92=sqrt(sum(ss.norm*(b18.g92[b.best[,1]]-mean.b18.g92)^2)/(w*sum(ss.norm)))

mean.b18.g39=sum(ss.norm*b18.g39[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.g39=sqrt(sum(ss.norm*(b18.g39[b.best[,1]]-mean.b18.g39)^2)/(w*sum(ss.norm)))

mean.b18.ps100198=sum(ss.norm*b18.ps100198[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.ps100198=sqrt(sum(ss.norm*(b18.ps100198[b.best[,1]]-mean.b18.ps100198)^2)/(w*sum(ss.norm)))


mean.b18.kapam=sum(ss.norm*b18.kapam[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.kapam=sqrt(sum(ss.norm*(b18.kapam[b.best[,1]]-mean.b18.kapam)^2)/(w*sum(ss.norm)))


mean.b18.storoen=sum(ss.norm*b18.storoen[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.storoen=sqrt(sum(ss.norm*(b18.storoen[b.best[,1]]-mean.b18.storoen)^2)/(w*sum(ss.norm)))

mean.b18.bourb=sum(ss.norm*b18.bourb[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.bourb=sqrt(sum(ss.norm*(b18.bourb[b.best[,1]]-mean.b18.bourb)^2)/(w*sum(ss.norm)))

mean.b18.blochn=sum(ss.norm*b18.blochn[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.blochn=sqrt(sum(ss.norm*(b18.blochn[b.best[,1]]-mean.b18.blochn)^2)/(w*sum(ss.norm)))


mean.b18.sondrem=sum(ss.norm*b18.sondrem[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.sondrem=sqrt(sum(ss.norm*(b18.sondrem[b.best[,1]]-mean.b18.sondrem)^2)/(w*sum(ss.norm)))

mean.b18.zi=sum(ss.norm*b18.zi[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.zi=sqrt(sum(ss.norm*(b18.zi[b.best[,1]]-mean.b18.zi)^2)/(w*sum(ss.norm)))

mean.b18.lambland1=sum(ss.norm*b18.lambland1[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.lambland1=sqrt(sum(ss.norm*(b18.lambland1[b.best[,1]]-mean.b18.lambland1)^2)/(w*sum(ss.norm)))

mean.b18.lambland2=sum(ss.norm*b18.lambland2[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.lambland2=sqrt(sum(ss.norm*(b18.lambland2[b.best[,1]]-mean.b18.lambland2)^2)/(w*sum(ss.norm)))

mean.b18.blaso=sum(ss.norm*b18.blaso[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b18.blaso=sqrt(sum(ss.norm*(b18.blaso[b.best[,1]]-mean.b18.blaso)^2)/(w*sum(ss.norm)))


###############################################################################
# B20
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)


mis.ps100=read.table(paste(work.fldr,"/rmse-ps100.txt",sep=""), skip=1)
mis.g92=read.table(paste(work.fldr,"/rmse-g92.txt",sep=""), skip=1)
mis.g39=read.table(paste(work.fldr,"/rmse-g39.txt",sep=""), skip=1)
mis.ps100198=read.table(paste(work.fldr,"/rmse-ps100198.txt",sep=""), skip=1)

mis.ps100=as.vector(mis.ps100[,1])
mis.g92=as.vector(mis.g92[,1])
mis.g39=as.vector(mis.g39[,1])
mis.ps100198=as.vector(mis.ps100198[,1])

mis.kapam=read.table(paste(work.fldr,"/rmse-kapam.txt",sep=""), skip=1)
mis.storoen=read.table(paste(work.fldr,"/rmse-storoen.txt",sep=""), skip=1)
mis.bourb=read.table(paste(work.fldr,"/rmse-bourb.txt",sep=""), skip=1)
mis.blochn=read.table(paste(work.fldr,"/rmse-blochn.txt",sep=""), skip=1)
mis.sondrem=read.table(paste(work.fldr,"/rmse-sondrem.txt",sep=""), skip=1)
mis.zi=read.table(paste(work.fldr,"/rmse-zi.txt",sep=""), skip=1)
mis.lambland1=read.table(paste(work.fldr,"/rmse-lambland1.txt",sep=""), skip=1)
mis.lambland2=read.table(paste(work.fldr,"/rmse-lambland2.txt",sep=""), skip=1)
mis.blaso=read.table(paste(work.fldr,"/rmse-blaso.txt",sep=""), skip=1)

mis.kapam=as.vector(mis.kapam[,1])
mis.storoen=as.vector(mis.storoen[,1])
mis.bourb=as.vector(mis.bourb[,1])
mis.blochn=as.vector(mis.blochn[,1])
mis.sondrem=as.vector(mis.sondrem[,1])
mis.zi=as.vector(mis.zi[,1])
mis.lambland1=as.vector(mis.lambland1[,1])
mis.lambland2=as.vector(mis.lambland2[,1])
mis.blaso=as.vector(mis.blaso[,1])


b20.ps100=mis.ps100
b20.g92=mis.g92
b20.g39=mis.g39
b20.ps100198=mis.ps100198
b20.kapam=mis.kapam
b20.storoen=mis.storoen
b20.bourb=mis.bourb
b20.blochn=mis.blochn
b20.sondrem=mis.sondrem
b20.zi=mis.zi
b20.lambland1=mis.lambland1
b20.lambland2=mis.lambland2
b20.blaso=mis.blaso


###############################################
# weighted mean
ss.norm=b.best[,17]
# median(b20.g92[b.best[,1]]

mean.b20.ps100=sum(ss.norm*b20.ps100[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.ps100=sqrt(sum(ss.norm*(b20.ps100[b.best[,1]]-mean.b20.ps100)^2)/(w*sum(ss.norm)))

mean.b20.g92=sum(ss.norm*b20.g92[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.g92=sqrt(sum(ss.norm*(b20.g92[b.best[,1]]-mean.b20.g92)^2)/(w*sum(ss.norm)))

mean.b20.g39=sum(ss.norm*b20.g39[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.g39=sqrt(sum(ss.norm*(b20.g39[b.best[,1]]-mean.b20.g39)^2)/(w*sum(ss.norm)))

mean.b20.ps100198=sum(ss.norm*b20.ps100198[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.ps100198=sqrt(sum(ss.norm*(b20.ps100198[b.best[,1]]-mean.b20.ps100198)^2)/(w*sum(ss.norm)))


mean.b20.kapam=sum(ss.norm*b20.kapam[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.kapam=sqrt(sum(ss.norm*(b20.kapam[b.best[,1]]-mean.b20.kapam)^2)/(w*sum(ss.norm)))


mean.b20.storoen=sum(ss.norm*b20.storoen[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.storoen=sqrt(sum(ss.norm*(b20.storoen[b.best[,1]]-mean.b20.storoen)^2)/(w*sum(ss.norm)))

mean.b20.bourb=sum(ss.norm*b20.bourb[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.bourb=sqrt(sum(ss.norm*(b20.bourb[b.best[,1]]-mean.b20.bourb)^2)/(w*sum(ss.norm)))

mean.b20.blochn=sum(ss.norm*b20.blochn[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.blochn=sqrt(sum(ss.norm*(b20.blochn[b.best[,1]]-mean.b20.blochn)^2)/(w*sum(ss.norm)))


mean.b20.sondrem=sum(ss.norm*b20.sondrem[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.sondrem=sqrt(sum(ss.norm*(b20.sondrem[b.best[,1]]-mean.b20.sondrem)^2)/(w*sum(ss.norm)))

mean.b20.zi=sum(ss.norm*b20.zi[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.zi=sqrt(sum(ss.norm*(b20.zi[b.best[,1]]-mean.b20.zi)^2)/(w*sum(ss.norm)))

mean.b20.lambland1=sum(ss.norm*b20.lambland1[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.lambland1=sqrt(sum(ss.norm*(b20.lambland1[b.best[,1]]-mean.b20.lambland1)^2)/(w*sum(ss.norm)))

mean.b20.lambland2=sum(ss.norm*b20.lambland2[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.lambland2=sqrt(sum(ss.norm*(b20.lambland2[b.best[,1]]-mean.b20.lambland2)^2)/(w*sum(ss.norm)))

mean.b20.blaso=sum(ss.norm*b20.blaso[b.best[,1]])/sum(ss.norm)
w=(length(ss.norm)-1)/length(ss.norm)
sd.b20.blaso=sqrt(sum(ss.norm*(b20.blaso[b.best[,1]]-mean.b20.blaso)^2)/(w*sum(ss.norm)))




####### plot
out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

plot.out=paste(work.fldr,"/Figures/fig4b-deglac-time.png", sep="")
png(plot.out, width=6, height=6, units="in", res=100)
#dev.new(width =6 , height = 6, units="in")
xlim=c(-15000,-8000)
ylim=c(-15000,0)
at.x=seq(-15000,-8000,by=1000)
label.x=seq(15,8,by=-1)
at.y=seq(-15000,0,by=1000)
label.y=seq(15,0,by=-1)

plot(xlim, ylim, type="n", ann=F, axes=T, xaxt="n", yaxt="n",cex.axis=1.5)
axis(2,at=at.y,labels=label.y,col="black",col.axis="black",las=1, cex.axis=1.3)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Observed deglaciation time (kyr ago)",side=1,line=2.75,cex=1.3)
mtext("Modelled deglaciation time (kyr ago)",side=2,line=2.75,cex=1.3)
abline(a=0, b=1, col="grey70")

# b18
points(ps100.age,mean.b18.ps100+ps100.age, bg="#390099",pch=24)
points(g92.age,mean.b18.g92+g92.age, bg="#390099",pch=24)
points(g39.age,mean.b18.g39+g39.age, bg="#390099", pch=24)
points(ps100198.age,mean.b18.ps100198+ps100198.age, bg="#390099",pch=24)

points(kapam.age,mean.b18.kapam+kapam.age, bg="#390099", pch=21)
points(storoen.age,mean.b18.storoen+storoen.age, bg="#390099", pch=21)
points(bourb.age,mean.b18.bourb+bourb.age, bg="#390099", pch=21)
points(blochn.age,mean.b18.blochn+blochn.age, bg="#390099", pch=21)
points(sondrem.age,mean.b18.sondrem+sondrem.age, bg="#390099",pch=21)
points(zi.age,mean.b18.zi+zi.age, bg="#390099", pch=21)
points(lambland1.age,mean.b18.lambland1+lambland1.age, bg="#390099", pch=21)
points(lambland2.age,mean.b18.lambland2+lambland2.age, bg="#390099", pch=21)
points(blaso.age,mean.b18.blaso+blaso.age, bg="#390099", pch=21)

# horiz error bars
arrows(x0=ps100.age-ps100.err, y0=mean.b18.ps100+ps100.age, x1=ps100.age+ps100.err, y1=mean.b18.ps100+ps100.age, code=3, length=0,lwd=1,lty=3, col=alpha.col("#390099",30))
arrows(x0=g92.age-g92.err, y0=mean.b18.g92+g92.age, x1=g92.age+g92.err, y1=mean.b18.g92+g92.age, code=3, length=0,lwd=1, lty=3, col=alpha.col("#390099",30))
arrows(x0=g39.age-g39.err, y0=mean.b18.g39+g39.age, x1=g39.age+g39.err, y1=mean.b18.g39+g39.age, code=3, length=0, lwd=1, lty=3, col=alpha.col("#390099",30))
arrows(x0=ps100198.age-ps100198.err, y0=mean.b18.ps100198+ps100198.age, x1=ps100198.age+ps100198.err, y1=mean.b18.ps100198+ps100198.age, code=3, length=0, lwd=1,lty=3, col=alpha.col("#390099",30))

arrows(x0=kapam.age-kapam.err, y0=mean.b18.kapam+kapam.age, x1=kapam.age+kapam.err, y1=mean.b18.kapam+kapam.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=storoen.age-storoen.err, y0=mean.b18.storoen+storoen.age, x1=storoen.age+storoen.err, y1=mean.b18.storoen+storoen.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=bourb.age-bourb.err, y0=mean.b18.bourb+bourb.age, x1=bourb.age+bourb.err, y1=mean.b18.bourb+bourb.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=blochn.age-blochn.err, y0=mean.b18.blochn+blochn.age, x1=blochn.age+blochn.err, y1=mean.b18.blochn+blochn.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=sondrem.age-sondrem.err, y0=mean.b18.sondrem+sondrem.age, x1=sondrem.age+sondrem.err, y1=mean.b18.sondrem+sondrem.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=zi.age-zi.err, y0=mean.b18.zi+zi.age, x1=zi.age+zi.err, y1=mean.b18.zi+zi.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=lambland1.age-lambland1.err, y0=mean.b18.lambland1+lambland1.age, x1=lambland1.age+lambland1.err, y1=mean.b18.lambland1+lambland1.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=lambland2.age-lambland2.err, y0=mean.b18.lambland2+lambland2.age, x1=lambland2.age+lambland2.err, y1=mean.b18.lambland2+lambland2.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=blaso.age-blaso.err, y0=mean.b18.blaso+blaso.age, x1=blaso.age+blaso.err, y1=mean.b18.blaso+blaso.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))

# vert error bars
arrows(x0=ps100.age, y0=mean.b18.ps100+ps100.age-sd.b18.ps100, x1=ps100.age, y1=mean.b18.ps100+ps100.age+sd.b18.ps100, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=g92.age, y0=mean.b18.g92+g92.age-sd.b18.g92, x1=g92.age, y1=mean.b18.g92+g92.age+sd.b18.g92, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=g39.age, y0=mean.b18.g39+g39.age-sd.b18.g39, x1=g39.age, y1=mean.b18.g39+g39.age+sd.b18.g39, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=ps100198.age, y0=mean.b18.ps100198+ps100198.age-sd.b18.ps100198, x1=ps100198.age, y1=mean.b18.ps100198+ps100198.age+sd.b18.ps100198, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))

arrows(x0=kapam.age, y0=mean.b18.kapam+kapam.age-sd.b18.kapam, x1=kapam.age, y1=mean.b18.kapam+kapam.age+sd.b18.kapam, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=storoen.age, y0=mean.b18.storoen+storoen.age-sd.b18.storoen, x1=storoen.age, y1=mean.b18.storoen+storoen.age+sd.b18.storoen, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=bourb.age, y0=mean.b18.bourb+bourb.age-sd.b18.bourb, x1=bourb.age, y1=mean.b18.bourb+bourb.age+sd.b18.bourb, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=blochn.age, y0=mean.b18.blochn+blochn.age-sd.b18.blochn, x1=blochn.age, y1=mean.b18.blochn+blochn.age+sd.b18.blochn, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=sondrem.age, y0=mean.b18.sondrem+sondrem.age-sd.b18.sondrem, x1=sondrem.age, y1=mean.b18.sondrem+sondrem.age+sd.b18.sondrem, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=zi.age, y0=mean.b18.zi+zi.age-sd.b18.zi, x1=zi.age, y1=mean.b18.zi+zi.age+sd.b18.zi, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=lambland1.age, y0=mean.b18.lambland1+lambland1.age-sd.b18.lambland1, x1=lambland1.age, y1=mean.b18.lambland1+lambland1.age+sd.b18.lambland1, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
arrows(x0=lambland2.age, y0=mean.b18.lambland2+lambland2.age-sd.b18.lambland2, x1=lambland2.age, y1=mean.b18.lambland2+lambland2.age+sd.b18.lambland2, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#390099",30))
#arrows(x0=blaso.age, y0=mean.b18.blaso+blaso.age-sd.b18.blaso, x1=blaso.age, y1=mean.b18.blaso+blaso.age+sd.b18.blaso, code=3, angle=90, length=0.05, lwd=1)


# b20
points(ps100.age,mean.b20.ps100+ps100.age, bg="#ffbd00", pch=24)
points(g92.age,mean.b20.g92+g92.age, bg="#ffbd00", pch=24)
points(g39.age,mean.b20.g39+g39.age, bg="#ffbd00", pch=24)
points(ps100198.age,mean.b20.ps100198+ps100198.age, bg="#ffbd00", pch=24)

points(kapam.age,mean.b20.kapam+kapam.age, bg="#ffbd00",pch=21)
points(storoen.age,mean.b20.storoen+storoen.age, bg="#ffbd00",pch=21)
points(bourb.age,mean.b20.bourb+bourb.age, bg="#ffbd00", pch=21)
points(blochn.age,mean.b20.blochn+blochn.age, bg="#ffbd00", pch=21)
points(sondrem.age,mean.b20.sondrem+sondrem.age, bg="#ffbd00", pch=21)
points(zi.age,mean.b20.zi+zi.age, bg="#ffbd00", pch=21)
points(lambland1.age,mean.b20.lambland1+lambland1.age, bg="#ffbd00",pch=21)
points(lambland2.age,mean.b20.lambland2+lambland2.age, bg="#ffbd00", pch=21)
points(blaso.age,mean.b20.blaso+blaso.age, bg="#ffbd00", pch=21)


# horiz error bars
arrows(x0=ps100.age-ps100.err, y0=mean.b20.ps100+ps100.age, x1=ps100.age+ps100.err, y1=mean.b20.ps100+ps100.age, code=3, length=0,lwd=1,lty=3, col=alpha.col("#ffbd00",30))
arrows(x0=g92.age-g92.err, y0=mean.b20.g92+g92.age, x1=g92.age+g92.err, y1=mean.b20.g92+g92.age, code=3, length=0,lwd=1, lty=3, col=alpha.col("#ffbd00",30))
arrows(x0=g39.age-g39.err, y0=mean.b20.g39+g39.age, x1=g39.age+g39.err, y1=mean.b20.g39+g39.age, code=3, length=0, lwd=1, lty=3, col=alpha.col("#ffbd00",30))
arrows(x0=ps100198.age-ps100198.err, y0=mean.b20.ps100198+ps100198.age, x1=ps100198.age+ps100198.err, y1=mean.b20.ps100198+ps100198.age, code=3, length=0, lwd=1,lty=3, col=alpha.col("#ffbd00",30))

arrows(x0=kapam.age-kapam.err, y0=mean.b20.kapam+kapam.age, x1=kapam.age+kapam.err, y1=mean.b20.kapam+kapam.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=storoen.age-storoen.err, y0=mean.b20.storoen+storoen.age, x1=storoen.age+storoen.err, y1=mean.b20.storoen+storoen.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=bourb.age-bourb.err, y0=mean.b20.bourb+bourb.age, x1=bourb.age+bourb.err, y1=mean.b20.bourb+bourb.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=blochn.age-blochn.err, y0=mean.b20.blochn+blochn.age, x1=blochn.age+blochn.err, y1=mean.b20.blochn+blochn.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=sondrem.age-sondrem.err, y0=mean.b20.sondrem+sondrem.age, x1=sondrem.age+sondrem.err, y1=mean.b20.sondrem+sondrem.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=zi.age-zi.err, y0=mean.b20.zi+zi.age, x1=zi.age+zi.err, y1=mean.b20.zi+zi.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=lambland1.age-lambland1.err, y0=mean.b20.lambland1+lambland1.age, x1=lambland1.age+lambland1.err, y1=mean.b20.lambland1+lambland1.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=lambland2.age-lambland2.err, y0=mean.b20.lambland2+lambland2.age, x1=lambland2.age+lambland2.err, y1=mean.b20.lambland2+lambland2.age, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=blaso.age-blaso.err, y0=mean.b20.blaso+blaso.age, x1=blaso.age+blaso.err, y1=mean.b20.blaso+blaso.age, code=3, angle=90, length=0.05, lwd=1,col=alpha.col("#ffbd00",30))

# vert error bars
arrows(x0=ps100.age, y0=mean.b20.ps100+ps100.age-sd.b20.ps100, x1=ps100.age, y1=mean.b20.ps100+ps100.age+sd.b20.ps100, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=g92.age, y0=mean.b20.g92+g92.age-sd.b20.g92, x1=g92.age, y1=mean.b20.g92+g92.age+sd.b20.g92, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=g39.age, y0=mean.b20.g39+g39.age-sd.b20.g39, x1=g39.age, y1=mean.b20.g39+g39.age+sd.b20.g39, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=ps100198.age, y0=mean.b20.ps100198+ps100198.age-sd.b20.ps100198, x1=ps100198.age, y1=mean.b20.ps100198+ps100198.age+sd.b20.ps100198, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))

arrows(x0=kapam.age, y0=mean.b20.kapam+kapam.age-sd.b20.kapam, x1=kapam.age, y1=mean.b20.kapam+kapam.age+sd.b20.kapam, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=storoen.age, y0=mean.b20.storoen+storoen.age-sd.b20.storoen, x1=storoen.age, y1=mean.b20.storoen+storoen.age+sd.b20.storoen, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=bourb.age, y0=mean.b20.bourb+bourb.age-sd.b20.bourb, x1=bourb.age, y1=mean.b20.bourb+bourb.age+sd.b20.bourb, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=blochn.age, y0=mean.b20.blochn+blochn.age-sd.b20.blochn, x1=blochn.age, y1=mean.b20.blochn+blochn.age+sd.b20.blochn, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=sondrem.age, y0=mean.b20.sondrem+sondrem.age-sd.b20.sondrem, x1=sondrem.age, y1=mean.b20.sondrem+sondrem.age+sd.b20.sondrem, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=zi.age, y0=mean.b20.zi+zi.age-sd.b20.zi, x1=zi.age, y1=mean.b20.zi+zi.age+sd.b20.zi, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=lambland1.age, y0=mean.b20.lambland1+lambland1.age-sd.b20.lambland1, x1=lambland1.age, y1=mean.b20.lambland1+lambland1.age+sd.b20.lambland1, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=lambland2.age, y0=mean.b20.lambland2+lambland2.age-sd.b20.lambland2, x1=lambland2.age, y1=mean.b20.lambland2+lambland2.age+sd.b20.lambland2, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))
arrows(x0=blaso.age, y0=mean.b20.blaso+blaso.age-sd.b20.blaso, x1=blaso.age, y1=mean.b20.blaso+blaso.age+sd.b20.blaso, code=3, angle=90, length=0.05, lwd=1, col=alpha.col("#ffbd00",30))

###################################
# core/morain names
text(ps100.age,-10200,"PS100/270", cex=0.6)
text(g92.age,-12600,"G92", cex=0.6)
text(g39.age,-14000,"G39", cex=0.6)
text(ps100198.age,-10500,"PS100/198", cex=0.6)

text(kapam.age,-11300,"Kap Amelie", cex=0.6)
text(storoen.age,-10500,"Storoen", cex=0.6)
text(bourb.age,-9000,"Bourbon Oer", cex=0.6)
text(blochn.age,-9200,"Bloch Nunatakker", cex=0.6)
text(sondrem.age,-2500,"Sondre Mellemland", cex=0.6)
text(zi.age,-6000,"ZI", cex=0.6)
text(lambland1.age,-4000,"Lambert Land1", cex=0.6)
text(lambland2.age,-7800,"Lambert Land2", cex=0.6)
text(blaso.age,-9800,"Blaso", cex=0.6)

dev.off()


### map

topo="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-4KM/GRL-4KM_TOPO-M17.nc"
htopo=raster(topo,varname="H_ice")
mask=raster(topo, varname="mask")
bed=raster(topo, varname="z_bed")
surf=raster(topo, varname="z_srf")
#nc=nc_open(topo)
lat2D=raster(topo,varname="lat2D")
lon2D=raster(topo,varname="lon2D")


#file="/home/titan/gwgk/gwgk005h/work/79N/data/topography/BedMachineGreenland-2021-04-20-NE.nc"
#bed=raster(file, varname="bed")
#thick=raster(file, varname="thickness")
#surf=raster(file, varname="surface")
#mask=raster(file, varname="mask")


# reproj bedmachine over yelmo
#bed.yelmo=projectRaster(bed, crs="+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )
#extent(bed.yelmo)=extent(lon2D)
#
#thick.yelmo=projectRaster(thick, crs="+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )
#extent(thick.yelmo)=extent(lon2D)
#
#surf.yelmo=projectRaster(surf, crs="+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )
#extent(surf.yelmo)=extent(lon2D)
#
#mask.yelmo=projectRaster(mask, crs="+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" )
#extent(mask.yelmo)=extent(lon2D)
#


# reproj to create lat lon
#bed.latlon=projectRaster(bed, crs="+init=EPSG:4326")
#dim(bed.latlon)=dim(bed)
#bed.latlon[]=bed[]
#
#g=as(bed.latlon, 'SpatialGridDataFrame')
#xcoor=coordinates(g)[,1]
#ycoor=coordinates(g)[,2]
#
## create a 2D map for x and y, with same extension as data
#lon2D=bed
#lat2D=bed
#lon2D[]=lat2D[]=NA
#
#values(lon2D)=xcoor
#values(lat2D)=ycoor

library(rgdal)

# core PS100/270-1 (https://doi.org/10.1594/PANGAEA.921185)
lat.ps100=79.497170
lon.ps100=-18.139670

# core DA17-NG-ST-08-092G (https://doi.pangaea.de/10.1594/PANGAEA.945987)
lat.g92=78.500900
lon.g92=-17.278517

# core DA17-NG-ST03-039G (https://doi.pangaea.de/10.1594/PANGAEA.943880) (Hansen et al., 2022, https://doi.org/10.1016/j.quascirev.2022.107704)
lat.g39=80.037167
lon.g39=-8.923183

# Core PS10/198 (Lloyd et al., 2023)
lat.ps100198=79.191233
lon.ps100198=-17.107167

#latlon.proj <- "+init=epsg:4326"
#bedma.proj <- "+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +rf=298.27940504282 +units=m +no_defs"
#yelmo.proj <- "+proj=stere +lat_0=90 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#
cores.latlon=data.frame(x=c(lon.ps100,lon.g92,lon.g39,lon.ps100198), y=c(lat.ps100,lat.g92,lat.g39,lat.ps100198))
#coordinates(cores.latlon) <- ~x+y
#proj4string(cores.latlon) <- CRS(latlon.proj)
#cores.bedma = spTransform(cores.latlon, CRS(bedma.proj))
cores.yelmo=data.frame(x=c(ps100.x,g92.x,g39.x,ps100198.x), y=c(ps100.y,g92.y,g39.y,ps100198.y))

# boulder/morains data

# OUTER COAST SITES 
# Bourbon Oer (Larsen et al., 2018)
lat.bourb=78.620
lon.bourb=-18.402

# Storoen (Larsen et al., 2018)
lat.storoen=78.066
lon.storoen=-19.091

# Kap Amelie (Larsen et al., 2018)
lat.kapam=77.544
lon.kapam=-19.131

# INNER (PRESENT-DAY) SITES
# Blaso (Larsen et al., 2018)
lat.blaso=79.636
lon.blaso=-23.104

# Lambert Land 1 (Larsen et al., 2018)
lat.lambland1=79.143
lon.lambland1=-21.41

# Lambert Land 2 (Larsen et al., 2018)
lat.lambland2= 79.1
lon.lambland2=-20.9

# Zachariae Istrom (Larsen et al., 2018)
lat.zi=78.782
lon.zi=-20.777

# Sondre Mellemland (Larsen et al., 2018)
lat.sondrem=78.067
lon.sondrem=-21.505

# Bloch Nunatakker (Larsen et al., 2018)
lat.blochn=79.6
lon.blochn=-19.524


morain.latlon=data.frame(x=c(lon.bourb,lon.storoen,lon.kapam,lon.blaso,lon.lambland1,lon.lambland2,lon.zi,lon.sondrem,lon.blochn), y=c(lat.bourb,lat.storoen,lat.kapam,lat.blaso,lat.lambland1,lat.lambland2,lat.zi,lat.sondrem,lat.blochn))
#coordinates(morain.latlon) <- ~x+y
#proj4string(morain.latlon) <- CRS(latlon.proj)
#morain.bedma = spTransform(morain.latlon, CRS(bedma.proj))
morain.yelmo=data.frame(x=c(bourb.x,storoen.x,kapam.x,blaso.x,lambland1.x,lambland2.x,zi.x,sondrem.x,blochn.x), y=c(bourb.y,storoen.y,kapam.y,blaso.y,lambland1.y,lambland2.y,zi.y,sondrem.y,blochn.y))

plot.out=paste(work.fldr,"/Figures/fig4a-map-morain-cores.png", sep="")
png(plot.out, width=4.5, height=4.5, units="in", res=100)
#dev.new(width =4.5 , height = 4.5, units="in")
#par(mfrow=c(1,2))
par(mar=c(1,1,1,0), oma=c(0,0,0,0))
surf[surf<1]=NA
col.bed=colorRampPalette(c("#ede0d4", "#e6ccb2", "#ddb892", "#b08968", "#7f5539", "#9c6644"))(100)
mask[mask<1]=NA
htopo[htopo<1]=NA
bed[bed < -400]= -400
topo=merge(surf,bed)
topo[topo>400]=400

#xlim=c(3e5,8e5)
#ylim=c(-1400000, -663875)
xlim=c(200,800)
ylim=c(-1400,-700)
plot(xlim,ylim, type="n", ann=F, axes=F, xaxt="n",yaxt="n", cex.axis=1.5)
plot(topo, add=T,col=alpha.col(col.bed,70),axes=F,legend=T)
#plot(surf, add=T, col=alpha.col(colorRampPalette(brewer.pal(9,"Blues"))(100),100), box=F, axes=F,legend=F)
#plot(surf, add=T, col=c("grey"),box=F, axes=F,legend=F)
plot(htopo, add=T, col="grey90", axes=F,legend=F)
contour(lat2D, nlevels=8, drawlabels=F, col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
contour(lon2D, nlevels=8, drawlabels=F,col=alpha.col("grey20", 40), lwd=0.5,lty=2,add=T)
#contour(mask, add=T, nlevels=1, drawlabels=F, lwd=2)
#points(coordinates(cores.bedma)[,1], coordinates(cores.bedma)[,2], pch=24, cex=1.5, bg="green")
#points(coordinates(morain.bedma)[,1], coordinates(morain.bedma)[,2], pch=21, cex=1.5, bg="green")
points(coordinates(cores.yelmo)[,1], coordinates(cores.yelmo)[,2], pch=24, cex=1.5, bg="green")
#points(coordinates(cores.bedma)[,1], coordinates(cores.bedma)[,2], pch=24, cex=1.5, bg="green")
points(coordinates(morain.yelmo)[,1], coordinates(morain.yelmo)[,2], pch=21, cex=1.5, bg="green")


text(coordinates(cores.yelmo)[1,1]+50, coordinates(cores.yelmo)[1,2]+20, "PS100/270", cex=0.8, font=2)
text(coordinates(cores.yelmo)[2,1]+30, coordinates(cores.yelmo)[2,2]-10, "G92", cex=0.8, font=2)
text(coordinates(cores.yelmo)[3,1], coordinates(cores.yelmo)[3,2]+30, "G39", cex=0.8, font=2)
text(coordinates(cores.yelmo)[4,1]+60, coordinates(cores.yelmo)[4,2]+10, "PS100/198", cex=0.8, font=2)
text(coordinates(morain.yelmo)[1,1]+50, coordinates(morain.yelmo)[1,2]+20, "Bourbon Oer", cex=0.8, font=2)
text(coordinates(morain.yelmo)[2,1]+50, coordinates(morain.yelmo)[2,2], "Storoen", cex=0.8, font=2)
text(coordinates(morain.yelmo)[3,1], coordinates(morain.yelmo)[3,2]-30, "Kap Amelie", cex=0.8, font=2)
text(coordinates(morain.yelmo)[4,1]-40, coordinates(morain.yelmo)[4,2], "Blaso", cex=0.8, font=2)
text(coordinates(morain.yelmo)[5,1]-80, coordinates(morain.yelmo)[5,2], "Lambert Land1", cex=0.8, font=2)
text(coordinates(morain.yelmo)[6,1]+60, coordinates(morain.yelmo)[6,2]+25, "Lambert Land2", cex=0.8, font=2)
text(coordinates(morain.yelmo)[7,1]-20, coordinates(morain.yelmo)[7,2]-15, "ZI", cex=0.8, font=2)
text(coordinates(morain.yelmo)[8,1], coordinates(morain.yelmo)[8,2]-20, "Sondre Mellemland", cex=0.8, font=2)
text(coordinates(morain.yelmo)[9,1]-90, coordinates(morain.yelmo)[9,2]+25, "Bloch Nunatakker", cex=0.8, font=2)

#text(coordinates(cores.bedma)[1,1]+20000, coordinates(cores.bedma)[1,2]+20000, "PS100/270", cex=0.5, font=2)
#text(coordinates(cores.bedma)[2,1]+15000, coordinates(cores.bedma)[2,2]+15000, "G92", cex=0.5, font=2)
#text(coordinates(cores.bedma)[3,1], coordinates(cores.bedma)[3,2]+20000, "G39", cex=0.5, font=2)
#text(coordinates(cores.bedma)[4,1]+20000, coordinates(cores.bedma)[4,2]+20000, "PS100/198", cex=0.5, font=2)
#text(coordinates(morain.bedma)[1,1], coordinates(morain.bedma)[1,2]-15000, "Bourbon Oer", cex=0.5, font=2)
#text(coordinates(morain.bedma)[2,1]+35000, coordinates(morain.bedma)[2,2], "Storoen", cex=0.5, font=2)
#text(coordinates(morain.bedma)[3,1], coordinates(morain.bedma)[3,2]-20000, "Kap Amelie", cex=0.5, font=2)
#text(coordinates(morain.bedma)[4,1]-25000, coordinates(morain.bedma)[4,2], "Blaso", cex=0.5, font=2)
#text(coordinates(morain.bedma)[5,1]-50000, coordinates(morain.bedma)[5,2], "Lambert Land1", cex=0.5, font=2)
#text(coordinates(morain.bedma)[6,1]+40000, coordinates(morain.bedma)[6,2]+15000, "Lambert Land2", cex=0.5, font=2)
#text(coordinates(morain.bedma)[7,1]-15000, coordinates(morain.bedma)[7,2]-10000, "ZI", cex=0.5, font=2)
#text(coordinates(morain.bedma)[8,1], coordinates(morain.bedma)[8,2]-20000, "Sondre Mellemland", cex=0.5, font=2)
#text(coordinates(morain.bedma)[9,1]-60000, coordinates(morain.bedma)[9,2]+15000, "Bloch Nunatakker", cex=0.5, font=2)


dev.off()


