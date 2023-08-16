# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

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
col=rev(plasma(nr.best))
col=rev(colorRampPalette(c("#390641","#6a994e","#a7c957","#bc4749"))(nr.best))
col=rev(mako(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")

#######################################################################################################################
# thin rates vs tau (first slope)
plot.out=paste(work.fldr,"/Figures/figS02-thin-rates-vs-taub-slopes.png", sep="")
png(plot.out, width = 10, height = 5, units = "in", res=100)
#dev.new(width = 10, height = 5, units="in")
par(mfrow=c(1,2))
xlim=c(40000,100000)
ylim=c(-50,0)
plot(xlim, ylim, type="n", xlab="Tb NEGIS (Pa)", ylab="NGRIP thinning rate (m/kyr)",main="Slope max-6 kyr")

print("first slope")
b=100
a=50
for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  min.i = which(time==-6000)

  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa 

  # plot
  points(taub.negis, thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
  nc_close(nc)
}

# Pearson correlation
print("pearson corr")
thin.rate=taub.negis=c()
for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
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
  max.i=which(zs.ngrip==max(zs.ngrip))
  min.i=which(time==-6000)

  thin.rate[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa 
  print(i)
}

res=cor.test(taub.negis, thin.rate, method="pearson")
text(55000,-10,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(55000,-15,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ taub.negis)
taub=seq(40000,100000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

# Thinning rates vs tau (second slope)
print("second slope")
xlim=c(40000,100000)
ylim=c(-20,10)
plot(xlim, ylim, type="n", xlab="Tb NEGIS (Pa)", ylab="NGRIP thinning rate (m/kyr)",main="Slope 6-0 kyr")

for(i in nr.best:1){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="") 
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  #max.i=which(zs.ngrip==max(zs.ngrip))
  #pd.i=length(time)
  max.i = which(time==-6000)
  min.i=which(time==0)

  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # calc beta averaged over NEGIS
  taub.hol=taub[,,max.i:min.i]
  taub.negis=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa 

  # plot
  points(taub.negis, thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
  nc_close(nc)
}

# Pearson correlation
print("pearson corr")
thin.rate=taub.negis=c()
for(i in 1:nr.best){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
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
  max.i=which(time==-6000)
  min.i=which(time==0)

  thin.rate[i]=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
  taub.hol=taub[,,max.i:min.i]
  taub.negis[i]=mean(taub.hol[basins>8.5 & basins<9.5])  # Pa 
  print(i)
}

res=cor.test(taub.negis, thin.rate, method="pearson")
text(70000,-15,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(70000,-18,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ taub.negis)
taub=seq(30000,80000)
lines(taub, taub*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()



