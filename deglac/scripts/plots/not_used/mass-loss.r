library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
col=rev(mako(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

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


#############################################################################################
# Thin rate vs mass Thin rate vs mass loss rate  

thin.rate=c()
mr=taub.sec=array(0, dim=c(nr.best,6))
for(i in 1:nr.best){
 i#file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")  
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  #max.i=which(dzs.ngrip==max(dzs.ngrip))
  max.i=38
  pd.i=length(time)
  min.i=pd.i
  thin.rate[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  #calc taub
  taub.ne=mean(taub[,,max.i:min.i][basins>1.5 & basins<2.5])  # Pa 
  taub.negis=mean(taub[,,max.i:min.i][basins>8.5 & basins<9.5])
  taub.n=mean(taub[,,max.i:min.i][basins>0.5 & basins<1.5])
  taub.nw=mean(taub[,,max.i:min.i][basins>7.5 & basins<8.5])
  taub.w=mean(taub[,,max.i:min.i][basins>6.5 & basins<7.5])
  taub.e=mean(taub[,,max.i:min.i][basins>2.5 & basins<3.5])
  taub.sw=mean(taub[,,max.i:min.i][basins>5.5 & basins<6.5])
  taub.sec[i,]=c(taub.e,taub.ne,taub.n,taub.nw,taub.w,taub.sw)

  # calc ice mass in Gt
  m.ne=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.ne=h.t[(basins>1.5 & basins<2.5) | (basins>8.5 & basins<9.5) ]
   h.ne=h.ne/1000  # km
   v.ne=sum(h.ne*8*8)  # km*km*km = km^3
   m.ne[t]=v.ne*0.9167  # km^3*Gt/km^3 = Gt
  }
  
  # calc ice mass in Gt
  m.n=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.n=h.t[basins>0.5 & basins<1.5]
   h.n=h.n/1000  # km
   v.n=sum(h.n*8*8)  # km*km*km = km^3
   m.n[t]=v.n*0.9167  # km^3*Gt/km^3 = Gt
  }

  # calc ice mass in Gt
  m.nw=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.nw=h.t[basins>7.5 & basins<8.5]
   h.nw=h.nw/1000  # km
   v.nw=sum(h.nw*8*8)  # km*km*km = km^3
   m.nw[t]=v.nw*0.9167  # km^3*Gt/km^3 = Gt
  }

  # calc ice mass in Gt
  m.w=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.w=h.t[basins>6.5 & basins<7.5]
   h.w=h.w/1000  # km
   v.w=sum(h.w*8*8)  # km*km*km = km^3
   m.w[t]=v.w*0.9167  # km^3*Gt/km^3 = Gt
  }

  # calc ice mass in Gt
  m.e=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.e=h.t[basins>2.5 & basins<3.5]
   h.e=h.e/1000  # km
   v.e=sum(h.e*8*8)  # km*km*km = km^3
   m.e[t]=v.e*0.9167  # km^3*Gt/km^3 = Gt
  }

  # calc ice mass in Gt
  m.sw=c()
  for(t in max.i:length(time)){
   h.t=h[,,t]
   h.sw=h.t[basins>5.5 & basins<6.5]
   h.sw=h.sw/1000  # km
   v.sw=sum(h.sw*8*8)  # km*km*km = km^3
   m.sw[t]=v.sw*0.9167  # km^3*Gt/km^3 = Gt
  }

  #mr.ne=(m.ne[max.i]-m.ne[min.i])/(time[max.i]-time[min.i])  # Gt/yr  
  #mr.e=(m.e[max.i]-m.e[min.i])/(time[max.i]-time[min.i])  
  #mr.n=(m.n[max.i]-m.n[min.i])/(time[max.i]-time[min.i])  
  #mr.nw=(m.nw[max.i]-m.nw[min.i])/(time[max.i]-time[min.i])  
  #mr.w=(m.w[max.i]-m.w[min.i])/(time[max.i]-time[min.i])  
  #mr.sw=(m.sw[max.i]-m.sw[min.i])/(time[max.i]-time[min.i])

  # mean rate for the last 9.5kyr
  max.i=which(time==-9600)
  min.i=length(time)
  mr.ne=(m.ne[max.i]-m.ne[min.i])/(time[max.i]-time[min.i])  # Gt/yr  
  mr.e=(m.e[max.i]-m.e[min.i])/(time[max.i]-time[min.i])
  mr.n=(m.n[max.i]-m.n[min.i])/(time[max.i]-time[min.i])
  mr.nw=(m.nw[max.i]-m.nw[min.i])/(time[max.i]-time[min.i])
  mr.w=(m.w[max.i]-m.w[min.i])/(time[max.i]-time[min.i])
  mr.sw=(m.sw[max.i]-m.sw[min.i])/(time[max.i]-time[min.i])

  mr[i,]=c(mr.e,mr.ne, mr.n, mr.nw, mr.w, mr.sw)
  #sector=c("E","NE","N","NW","W")

  print(i)
}

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-NE.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="NE mass rate (Gt/yr)")
points(rev(mr[,2]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,2], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,2])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-E.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="E mass rate (Gt/yr)")
points(rev(mr[,1]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,1], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,1])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-N.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
#xlim=c(0,1e7)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="N mass rate (Gt/yr)")
points(rev(mr[,3]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,3], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,3])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-NW.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
#xlim=c(0,1e7)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="NW mass rate (Gt/yr)")
points(rev(mr[,4]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,4], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,4])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-W.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
#xlim=c(0,1e7)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="W mass rate (Gt/yr)")
points(rev(mr[,5]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,5], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,5])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-mass-loss-SW.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(-30,10)
#xlim=c(0,1e7)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="SW mass rate (Gt/yr)")


# briner202 mass loss (Gt/century) (for last 10kyr)
data=read.table("/home/titan/gwgk/gwgk005h/work/ngrip/data/briner2020/mass-loss.txt", skip=144)
#data$V11=rowMeans(data[,c('V2','V3','V4','V5','V6','V7','V8','V9','V10')], na.rm=TRUE)
meanv2=mean(data$V2)
meanv3=mean(data$V3)
meanv4=mean(data$V4)
meanv5=mean(data$V5)
meanv6=mean(data$V6)
meanv7=mean(data$V7)
meanv8=mean(data$V8)
meanv9=mean(data$V9)
meanv10=mean(data$V10)
mean.mass.rate=c(meanv2,meanv3,meanv4,meanv5,meanv6,meanv7,meanv8,meanv9,meanv10)/100  # mass loss/year
briner.max=min(mean.mass.rate)
briner.min=max(mean.mass.rate)
#briner.max=-4.78
#briner.min=3.715

rect(briner.max,-20,briner.min,10,col="pink", border=NA)
points(rev(mr[,6]),rev(thin.rate), col="black",bg=rev(mycol), pch=21, cex=1)
abline(v=0, lty=2, col="grey50")
# correlation
res=cor.test(mr[,6], thin.rate, method="pearson")
text(-20,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-20,-25,paste("p = ",signif(res$p.value, digits=4),sep=""))
# linear regression
lin.reg=lm(thin.rate ~ mr[,6])
mass=seq(-80,80)
lines(mass, mass*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

#############################################################################
# boxplot

data=as.data.frame(mr)
colnames(data)=c("E","NE","N","NW","W","SW")

plot.out=paste(work.fldr,"/Figures/massloss-boxplot.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
boxplot(data, col=c("green","green","grey80","grey80","green","grey80"), xlab="GrIS Sector", ylab="Mass change rate (Gt/yr)")
grid()
legend("bottomright", legend=c("Corr. with thinning","no correlation"),col=c("green","grey"),bty="n",pch=20,pt.cex=3,cex=1,horiz=F)
dev.off()

