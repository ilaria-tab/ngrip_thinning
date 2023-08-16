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
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs (ngrip max thinning)
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,19],breaks = nr.best))] # mycol is wrt to thinning

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

###################################################################################################################

# Thin rates vs taub (NEGIS)

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))

## Pearson correlation
thin.best=taub.negis.best=thin=taub.negis=c()
for(i in 1:nr){
#   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b[i,1]-1,"/yelmo2D.nc",sep="")
 
  # does the file exist?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    thin[i]=0
    taub.negis[i]=60000  
    next
  }

  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  #time=time-1950
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip[15:nbands]==max(zs.ngrip[15:nbands]))
  pd.i=length(time)
  min.i=pd.i
  #min.i=62  

  # has the sim crashed?
  #nc = nc_open(file)
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    thin[i]=0
    taub.negis[i]=60000
    next
  }

  thin[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa  
  print(i)
}
## write results to a file#
file.out=paste(work.fldr,"/scripts/plots/thin.txt",sep="")
write.table(thin, file.out, quote = F, sep = " ", row.names=F)
file.out=paste(work.fldr,"/scripts/plots/taub-negis.txt",sep="")
write.table(taub.negis, file.out, quote = F, sep = " ", row.names=F)
#
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")  
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  #time=time-1950
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip[15:nbands]==max(zs.ngrip[15:nbands]))
  pd.i=length(time)
  min.i=pd.i
  #min.i = 62 # 5kyr ago  

  thin.best[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000 # m/kyr
  taub.hol.best=taub[,,max.i:min.i]
  taub.negis.best[i]=mean(taub.hol.best[basins>8.5 & basins<9.5])  # Pa  
  print(i)
}
### write results to a file
file.out=paste(work.fldr,"/scripts/plots/thin-best.txt",sep="")
write.table(thin.best, file.out, quote = F, sep = " ", row.names=F)
file.out=paste(work.fldr,"/scripts/plots/taub-negis-best.txt",sep="")
write.table(taub.negis.best, file.out, quote = F, sep = " ", row.names=F)




########################################################################################################

# READ files #

thin.best=read.table(paste(work.fldr,"/scripts/plots/thin-best.txt",sep=""), skip=1)
thin=read.table(paste(work.fldr,"/scripts/plots/thin.txt",sep=""), skip=1)
taub.negis.best=read.table(paste(work.fldr,"/scripts/plots/taub-negis-best.txt",sep=""), skip=1)
taub.negis=read.table(paste(work.fldr,"/scripts/plots/taub-negis.txt",sep=""), skip=1)

taub.negis=unlist(taub.negis, use.names=FALSE)
taub.negis.best=unlist(taub.negis.best, use.names=FALSE)
thin=unlist(thin, use.names=FALSE)
thin.best=unlist(thin.best, use.names=FALSE)

# calculate density of points
#library(MASS)
#z=kde2d(taub.negis, thin, n=50)

# which is the best simulation?
#best.sim=which(b.best[1,11]==b[,11])
best.sim=b.best[1,1]

# plot
plot.out=paste(work.fldr,"/Figures/fig5b-thin-rates-vs-taub-negis.png", sep="")
png(plot.out, width = 4, height = 5, units = "in",res=100)
#dev.new(width = 4, height = 5, units="in")
par(mar = c(0,0,0,0), oma = c(3.3,3.5,1,2.7),mgp=c(0,0.5,0))
#par(mar = c(0,0,0,0), oma = c(3.5,3.5,1,3))
#par(mar = c(4, 4, 1, 3) + 0.3)
at.x=seq(40000,90000,by=10000)
label.x=seq(40000,90000,by=10000)
at.y=seq(-20,5, by=5)
label.y=seq(-20,5,by=5)
xlim=c(40000,90000)
ylim=c(-40,5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n",cex.axis=1)

#contour(z, drawlabels=F, nlevels=7, col="grey75", add=T, lwd=1)
points(taub.negis, thin, pch=20, cex=0.3, col="grey30")
points(rev(taub.negis.best), rev(thin.best), pch=21, cex=1, col="black", bg=rev(mycol))
points(taub.negis[best.sim],thin[best.sim], pch=24, cex=1,col="blue", bg=mycol, lwd=2.5)

## thin rate Vinther 2009
max.v09=which(ds.ngrip.v09==max(ds.ngrip.v09))
pd.v09=length(t.ds.v09)
thin.rate.v09=(ds.ngrip.v09[max.v09]-ds.ngrip.v09[pd.v09])/(t.ds.v09[max.v09]-t.ds.v09[pd.v09])*1000
abline(h=thin.rate.v09, lwd=1,lty=3, col="black")

# linear regression BEST
res=cor.test(taub.negis.best, thin.best, method="pearson")
res.negis=res$estimate[[1]]
p.negis=res$p.value
text(80000,-10,paste("R = ",round(res$estimate[[1]], digits=2),sep=""),cex=1)
text(80000,-12,paste("p = ",signif(res$p.value, digits=2),sep=""), cex=0.8)
lin.reg=lm(thin.best ~ taub.negis.best)
taub=seq(10000,110000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30", lty=1)

axis(2,at=at.y,labels=label.y,col="black",col.axis="black",las=0, cex.axis=1,tck=-0.02)
mtext("					NGRIP thinning rate (m/kyr)",side=2,line=2.4,cex=1)
axis(1,at=at.x,labels=label.x,col="black",col.axis="black",las=1, cex.axis=1,tck=-0.015)
mtext(expression(paste("NEGIS ",tau[b], " (Pa)")),side=1,line=2.2,cex=1)

# plot skill score vs NEGIS beta (color is wrt ngrip thinning)
par(new=TRUE)
# taub vs score
plot(taub.negis, b[,17], pch=20, cex=0.3, col="grey30", xlim=c(40000,90000), ylim=c(0,2.8), axes=F, xlab="", ylab="")
# taub best vs score.best
points(taub.negis.best, b.best[,17], pch=21, cex=1, col="black", bg=mycol)
points(taub.negis[best.sim],b[best.sim,17], pch=24, cex=1,col="blue", bg=mycol, lwd=2.5)
#axis(side = 4, at = pretty(range(c(0,1))), cex.axis=1, tck=-0.02)
axis(side = 4, at = round(seq(0,1,by=0.5), digits=1), cex.axis=1, tck=-0.02)
mtext("Score                                              		", cex=1,side = 4, line = 1.5)

dev.off()

