# total score

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
source(paste(work.fldr,"/scripts/read_LHS.r",sep=""))

ngrip.max=read.table(paste(work.fldr,"/ngrip-max.txt",sep=""), skip=1)
ngrip.max=as.vector(ngrip.max[,1])

# new method with median of rmse for each misfit
mis.h=read.table(paste(work.fldr,"/rmse-h.txt",sep=""), skip=1)
#mis.s=read.table(paste(work.fldr,"/rmse-s.txt",sep=""), skip=1)
mis.cores=read.table(paste(work.fldr,"/rmse-cores.txt",sep=""), skip=1)
mis.h.ne=read.table(paste(work.fldr,"/rmse-h-ne.txt",sep=""), skip=1)
mis.u=read.table(paste(work.fldr,"/rmse-u.txt",sep=""), skip=1)
mis.ngrip=read.table(paste(work.fldr,"/rmse-ngrip.txt",sep=""), skip=1)
#mis.79.margin=read.table(paste(work.fldr,"/rmse-79-margin.txt",sep=""), skip=1)
mis.area.pd=read.table(paste(work.fldr,"/rmse-area-pd.txt",sep=""), skip=1)
mis.area.lgm=read.table(paste(work.fldr,"/rmse-area-lgm.txt",sep=""), skip=1)
mis.area.ne=read.table(paste(work.fldr,"/rmse-area-pd-ne.txt",sep=""), skip=1)
#mis.area.sw=read.table(paste(work.fldr,"/rmse-area-pd-sw.txt",sep=""), skip=1)
mis.u.ne=read.table(paste(work.fldr,"/rmse-u-ne.txt",sep=""), skip=1)
mis.ps100=read.table(paste(work.fldr,"/rmse-ps100.txt",sep=""), skip=1)
mis.g92=read.table(paste(work.fldr,"/rmse-g92.txt",sep=""), skip=1)
mis.g39=read.table(paste(work.fldr,"/rmse-g39.txt",sep=""), skip=1)
mis.ps100198=read.table(paste(work.fldr,"/rmse-ps100198.txt",sep=""), skip=1)
mis.marine=read.table(paste(work.fldr,"/rmse-marine.txt",sep=""), skip=1)
mis.inner.coast=read.table(paste(work.fldr,"/rmse-inner-coast.txt",sep=""), skip=1)
mis.outer.coast=read.table(paste(work.fldr,"/rmse-outer-coast.txt",sep=""), skip=1)
mis.ng=read.table(paste(work.fldr,"/rmse-ng.txt",sep=""), skip=1)

mis.h=as.vector(mis.h[,1])
mis.cores=as.vector(mis.cores[,1])
mis.h.ne=as.vector(mis.h.ne[,1])
mis.u=as.vector(mis.u[,1])
mis.ngrip=as.vector(mis.ngrip[,1])
#mis.79.margin=as.vector(mis.79.margin[,1])
mis.area.pd=as.vector(mis.area.pd[,1])
mis.area.lgm=as.vector(mis.area.lgm[,1])
mis.area.ne=as.vector(mis.area.ne[,1])
#mis.area.sw=as.vector(mis.area.sw[,1])
mis.u.ne=as.vector(mis.u.ne[,1])
mis.ps100=as.vector(mis.ps100[,1])
mis.g92=as.vector(mis.g92[,1])
mis.g39=as.vector(mis.g39[,1])
mis.ps100198=as.vector(mis.ps100198[,1])
mis.marine=as.vector(mis.marine[,1])
mis.inner.coast=as.vector(mis.inner.coast[,1])
mis.outer.coast=as.vector(mis.outer.coast[,1])
mis.ng=as.vector(mis.ng[,1])

s.h=exp(-mis.h/median(mis.h, na.rm=T))
s.cores=exp(-mis.cores/median(mis.cores, na.rm=T))
s.h.ne=exp(-mis.h.ne/median(mis.h.ne, na.rm=T))
s.u=exp(-mis.u/median(mis.u, na.rm=T))
s.ngrip=exp(-mis.ngrip/median(mis.ngrip, na.rm=T))
#s.79.margin=exp(-mis.79.margin/median(mis.79.margin, na.rm=T))
s.area.pd=exp(-mis.area.pd/median(mis.area.pd, na.rm=T))
s.area.lgm=exp(-mis.area.lgm/median(mis.area.lgm, na.rm=T))
s.area.ne=exp(-mis.area.ne/median(mis.area.ne, na.rm=T))
#s.area.sw=exp(-mis.area.sw/median(mis.area.sw, na.rm=T))
s.u.ne=exp(-mis.u.ne/median(mis.u.ne, na.rm=T))
#s.ps100=exp(-mis.ps100/median(mis.ps100, na.rm=T))
#s.g92=exp(-mis.g92/median(mis.g92, na.rm=T))
#s.g39=exp(-mis.g39/median(mis.g39, na.rm=T))
s.marine=exp(-mis.marine/median(mis.marine, na.rm=T))
s.inner.coast=exp(-mis.inner.coast/median(mis.inner.coast, na.rm=T))
s.outer.coast=exp(-mis.outer.coast/median(mis.outer.coast, na.rm=T))
s.ng=exp(-mis.ng/median(mis.ng, na.rm=T))
s.ng[s.ng>0.3]=1

ss = s.h*s.area.pd*s.area.lgm*s.ngrip*s.area.ne*s.cores*s.h.ne*s.u.ne*s.u*s.outer.coast*s.marine*s.ng*s.inner.coast
#ss = s.h * s.area.pd * s.ngrip * s.79.margin * s.area.ne * s.cores * s.h.ne * s.g92 * s.g39 *s.ps100 * s.u.ne *s.u
#ss = s.h * s.area.pd * s.ngrip * s.79.margin * s.cores  * s.g92 * s.g39 *s.ps100

# score normalization between 0 and 1
min.s=min(ss, na.rm=T)
max.s=max(ss, na.rm=T)
ss.norm=(ss-min.s)/(max.s-min.s)
ss.norm[is.na(ss.norm)==T]=0

b=cbind(a,ss.norm)
# best simulations in order
#b[order(b[,6]),]

#igood=which(mis.cores <20)
igood=which(mis.cores <50)
#igood=which(mis.cores < 100)
class=seq(1,nr)
class[igood]=1000
class[class != 1000]=0

# simulation number 
sim=seq(1,nr)
#b=cbind(sim,a,ss.norm,class,ngrip.max,s.u,s.h, s.area.pd,s.ngrip,s.79.margin,s.area.ne, s.cores,s.h.ne,s.u.ne)
b=cbind(sim,a,ss.norm,class,ngrip.max)
b.ord=b[order(b[,17], decreasing=T),]

best=c()
j=1
for(i in 1:nr){
 #if(b.ord[i,18]==1000 & b.ord[i,17] >= 0.001){
 if(b.ord[i,18]==1000 & b.ord[i,17] > 0){
 #if(b.ord[i,10] >= 0.5){ 
   best[j]=i
   j=j+1
 }else{
  # break
  j=j
  next
 }
}
b.best=b.ord[best,]
#b.best=b.best[1:34,]
nr.best=nrow(b.best)

itmb=round(as.vector(b.best[,2]),digits=4)
itmc=round(as.vector(b.best[,3]),digits=4)
q = round(as.vector(b.best[,4]),digits=4)
z0 = round(as.vector(b.best[,5]),digits=4)
neff = round(as.vector(b.best[,6]),digits=4)
enh = round(as.vector(b.best[,7]),digits=4)
k = round(as.vector(b.best[,8]),digits=4)
tau = round(as.vector(b.best[,9]),digits=4)
lith=round(as.vector(b.best[,10]),digits=4)
kt = round(as.vector(b.best[,11]),digits=4)
negis_c = round(as.vector(b.best[,12]),digits=4)
negis_s = round(as.vector(b.best[,13]),digits=4)
negis_n = round(as.vector(b.best[,14]),digits=4)
fp = round(as.vector(b.best[,15]),digits=4)
negis1=round(as.vector(b.best[,16]),digits=4)

itmb=sprintf('%.4f',itmb)
itmc=sprintf('%.4f',itmc)
q=sprintf('%.4f',q)
z0=sprintf('%.4f',z0)
neff=sprintf('%.4f',neff)
enh=sprintf('%.4f',enh)
k=sprintf('%.4f',k)
tau=sprintf('%.4f',tau)
lith=sprintf('%.4f',lith)
kt=sprintf('%.4f',kt)
negis_c=sprintf('%.4f',negis_c)
negis_s=sprintf('%.4f',negis_s)
negis_n=sprintf('%.4f',negis_n)
fp=sprintf('%.4f',fp)
negis1=sprintf('%.4f',negis1)

## discarded runs
#worst=c()
#j=1
#for(i in 1:nr){
#  if(b.ord[i,14]==0){
#  # if(b.ord[i,10] < 0.5){
#   worst[j]=i
#   j=j+1 
#  }
#}
#b.worst=b.ord[worst,]
#nr.worst=nrow(b.worst)
#
#itmb.worst=round(as.vector(b.worst[,2]),digits=4)
#itmc.worst=round(as.vector(b.worst[,3]),digits=4)
#q.worst = round(as.vector(b.worst[,4]),digits=4)
#z0.worst = round(as.vector(b.worst[,5]),digits=4)
#neff.worst = round(as.vector(b.worst[,6]),digits=4)
#enh.worst = round(as.vector(b.worst[,7]),digits=4)
#k.worst = round(as.vector(b.worst[,8]),digits=4)
#tau.worst = round(as.vector(b.worst[,9]),digits=4)
#kt.worst = round(as.vector(b.worst[,10]),digits=4)
#negis_c.worst=round(as.vector(b.worst[,11]),digits=4)
#negis_s.worst=round(as.vector(b.worst[,12]),digits=4)
#negis_n.worst=round(as.vector(b.worst[,13]),digits=4)
#
#itmb.worst=sprintf('%.4f',itmb.worst)
#itmc.worst=sprintf('%.4f',itmc.worst)
#q.worst=sprintf('%.4f',q.worst)
#z0.worst=sprintf('%.4f',z0.worst)
#neff.worst=sprintf('%.4f',neff.worst)
#enh.worst=sprintf('%.4f',enh.worst)
#k.worst=sprintf('%.4f',k.worst)
#tau.worst=sprintf('%.4f',tau.worst)
#kt.worst=sprintf('%.4f',kt.worst)
#negis_c.worst=sprintf('%.4f',negis_c.worst)
#negis_s.worst=sprintf('%.4f',negis_s.worst)
#negis_n.worst=sprintf('%.4f',negis_n.worst)
