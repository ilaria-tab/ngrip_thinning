# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90  # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
#col=rev(plasma(nr.best))
col=rev(mako(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")


##################################################################################################################

## this is a repetition of Fig3b (a)
## Thin vs taub (NEGIS)
#plot.out=paste(work.fldr,"/Figures/thin-vs-taub-negis.png", sep="")
#png(plot.out, width = 5, height = 5, units = "in", res=100)
##dev.new(width = 5, height = 5, units="in")
#xlim=c(30000,100000)
#ylim=c(-30,10)
#plot(xlim, ylim, type="n", xlab="Tb NEGIS (Pa)", ylab="NGRIP thinning (m)")
#
## Pearson correlation
#thin=taub.negis=c()
#for(i in 1:nr.best){
#   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
#  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
#  nc=nc_open(file)
#  basins=ncvar_get(nc,"basins")
#  beta=ncvar_get(nc,"beta")
#  taub=ncvar_get(nc,"taub")
#  u=ncvar_get(nc,"uxy_s")
#  zs=ncvar_get(nc,"z_srf")
#  time=ncvar_get(nc,"time")
 # nc_close(nc)
 # # calc thinning rate
 # zs.ngrip=zs[ngrip.i,ngrip.j,]
 # #max.i=which(zs.ngrip==max(zs.ngrip))
 # max.i=which(zs.ngrip[11:nbands]==max(zs.ngrip[11:nbands]))
 # pd.i=length(time)
 # min.i=pd.i
#
#  thin[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
#  #thin[i]=b.best[i,1]  # NGRIP thinning
#  taub.hol=taub[,,max.i:min.i]
#  taub.negis[i]=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa  
#}
#
#points(taub.negis, thin, pch=21, cex=1, col="black", bg=mycol)
#res=cor.test(taub.negis, thin, method="pearson")
#res.negis=res$estimate[[1]]
#p.negis=res$p.value
#text(40000,-15,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
#text(40000,-20,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
#
## linear regression
#lin.reg=lm(thin ~ taub.negis)
#taub=seq(20000,100000)
#lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
#dev.off()
#
##########################################################################################################################

# Thin rates vs taub (NORTH)
plot.out=paste(work.fldr,"/Figures/figS03a-thin-rates-vs-taub-north.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(30000,75000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", xlab="Tb NORTH (Pa)", ylab="Mean NGRIP thinning rate (m/kyr)")

thin=taub.negis=c()
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  #max.i=which(time==-10000)
  pd.i=length(time)
  min.i=pd.i
  thin[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000 # m/kyr

  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>0.5 & basins<1.5])  # Pa 
}

points(taub.negis, thin, col="black",bg=rev(mycol),pch=21, cex=1)
res=cor.test(taub.negis, thin, method="pearson")
res.n=res$estimate[[1]]
p.n=res$p.value
text(60000,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(60000,-25,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ taub.negis)
taub=seq(20000,75000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()

##########################################################################################################################

# Thin rates vs taub (NORTH-WEST)
plot.out=paste(work.fldr,"/Figures/figS03b-thin-rates-vs-taub-north-west.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(30000,75000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", xlab="Tb NORTH-WEST (Pa)", ylab="Mean NGRIP thinning rate (m/kyr)")

thin=taub.negis=c()
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>7.5 & basins<8.5])  # Pa 

  # plot
 # points(taub.negis, thin, col="black",bg=mycol[i],pch=21, cex=1)
}

points(taub.negis, thin, col="black",bg=rev(mycol),pch=21, cex=1)
res=cor.test(taub.negis, thin, method="pearson")
res.nw=res$estimate[[1]]
p.nw=res$p.value
text(60000,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(60000,-25,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ taub.negis)
taub=seq(20000,75000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()

########################################################################################################################
# Thin rates vs taub (WEST)
plot.out=paste(work.fldr,"/Figures/figS03c-thin-rates-vs-taub-west.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(50000,100000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", xlab="Tb WEST (Pa)", ylab="Mean NGRIP thinning rate (m/kyr)")

thin=taub.negis=c()
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>6.5 & basins<7.5])  # Pa 

  # plot
  #points(taub.negis, thin, col="black",bg=mycol[i],pch=21, cex=1)
}

points(taub.negis, thin, col="black",bg=rev(mycol),pch=21, cex=1)
res=cor.test(taub.negis, thin, method="pearson")
res.w=res$estimate[[1]]
p.w=res$p.value
text(60000,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(60000,-25,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ taub.negis)
taub.new=seq(50000,100000)
lines(taub.new, taub.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()

########################################################################################################################

# Thin rates vs taub (EAST)
plot.out=paste(work.fldr,"/Figures/figS03d-thin-rates-vs-taub-east.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(20000,75000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", xlab="Tb EAST (Pa)", ylab="Mean NGRIP thinning rate (m/kyr)")

thin=taub.negis=c()
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>2.5 & basins<3.5])  # Pa 
}

points(taub.negis, thin, col="black",bg=rev(mycol),pch=21, cex=1)
res=cor.test(taub.negis, thin, method="pearson")
res.e=res$estimate[[1]]
p.e=res$p.value
text(60000,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(60000,-25,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ taub.negis)
taub=seq(20000,75000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()

############################################################################################################
# Thin rates vs taub (NORTH-EAST)
plot.out=paste(work.fldr,"/Figures/figS03e-thin-rates-vs-taub-northeast.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(30000,75000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", xlab="Tb NORTH-EAST (Pa)", ylab="Mean NGRIP thinning rate (m/kyr)")

thin=taub.negis=c()
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  thin[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000 # m/kyr
  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub[(basins>1.5 & basins<2.5) | (basins>8.5 & basins<9.5)])  # Pa 
}

points(taub.negis, thin, col="black",bg=rev(mycol),pch=21, cex=1)
res=cor.test(taub.negis, thin, method="pearson")
res.ne=res$estimate[[1]]
text(60000,-20,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(60000,-25,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ taub.negis)
taub=seq(20000,75000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()


