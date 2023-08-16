library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,9],breaks = nr.best))]

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

############################################################################
# U NEGIS
plot.out=paste(work.fldr,"/Figures/thin-rates-vs-u-negis.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(50,200)
ylim=c(-20,0)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="U NEGIS (m/yr)")

for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
 
  # U for the Holocene
  u.hol=u[,,max.i:min.i]
  u.ne=mean(u.hol[basins==9])

  points(u.ne,thin.rate,col="black",bg=mycol[i],pch=21, cex=1.3)
}

# Pearson correlation
thin.rate=u.ne=c()
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  h = ncvar_get(nc,"H_ice")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i = which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i=pd.i

  thin.rate[i]= (zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
  u.hol=u[,,max.i:min.i]
  u.ne[i]=mean(u.hol[basins==9])
}

res=cor.test(u.ne, thin.rate, method="pearson")
text(70,-17,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(70,-18,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ u.ne)
u=seq(-40,200)
lines(u, u*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()
##########################################################################################################################
# U EGRIP
plot.out=paste(work.fldr,"/Figures/thin-rates-vs-u-egrip.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,100)
ylim=c(-20,0)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="U EGRIP (m/yr)")

for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
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
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  u.ne=mean(u[egrip.i, egrip.j,max.i:min.i])
  
  points(u.ne,thin.rate,col="black",bg=mycol[i],pch=21, cex=1.3)
}

# Pearson correlation
thin.rate=u.ne=c()
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  h = ncvar_get(nc,"H_ice")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)
  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i = which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i=pd.i

  thin.rate[i]= (zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
  u.ne[i]=mean(u[egrip.i, egrip.j,max.i:min.i])
}

res=cor.test(u.ne, thin.rate, method="pearson")
text(10,-17,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(10,-18,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ u.ne)
u=seq(-40,200)
lines(u, u*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

dev.off()
