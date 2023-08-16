library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(viridis)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
#col=rev(plasma(nr))
col=rev(mako(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]


# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")


############################################################################
# Plot
plot.out=paste(work.fldr,"/Figures/figS04-thin-rates-vs-all_param.png", sep="")
png(plot.out, width = 12, height = 9, units = "in", res=100)

#dev.new(width = 12, height = 9, units="in")
#par(mfrow=c(2,4), mar=c(5,5,1,1), cex.axis=1.2, cex.lab=1.2)
par(mfrow=c(3,4), mar=c(5,5,1,1), cex.axis=1.2, cex.lab=1.2)



############################################################################
# cf_NEGIS

xlim=c(0,1)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="cf_NEGIS")
#axis(1,cex.axis=1.2)
#axis(2, cex.axis=1.2)

for(i in nr.best:1){
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
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(negis_c[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}

# Pearson correlation
thin=c()
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
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)

  thin[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000 # m/kyr
}

res=cor.test(as.numeric(negis_c), thin, method="pearson")
text(0.8,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(0.8,-28,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(negis_c))
cf=seq(0,1)
lines(cf, cf*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

###########################################################################
# q
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,1)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="q")

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(q[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}

res=cor.test(as.numeric(q), thin, method="pearson")
text(0.1,-18,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(0.1,-22,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(q))
q.exp=seq(0,1)
lines(q.exp, q.exp*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

###########################################################################
# z0
#dev.new(width = 5, height = 5, units="in")
xlim=c(-1000,0)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="z0")

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

   # plot
  points(z0[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(z0), thin, method="pearson")
text(-400,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-400,-28,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(z0))
z0.new=seq(-800,100)
lines(z0.new, z0.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")

###########################################################################
# delta
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,0.06)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(delta))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(neff[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(neff), thin, method="pearson")
text(0.01,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(0.01,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(neff))
delta.new=seq(0,0.06,0.001)
lines(delta.new, delta.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")



###########################################################################
# itm_b
xlim=c(-3.5,-2)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="itm_b")
#axis(1,cex.axis=1.2)
#axis(2, cex.axis=1.2)

for(i in nr.best:1){
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
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr
  
  # plot
  points(itmb[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}

# Pearson correlation
thin=c()
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
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)

  thin[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000 # m/kyr
}

res=cor.test(as.numeric(itmb), thin, method="pearson")
text(-3,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-3,-28,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(itmb))
itmb.new=seq(-3.5,-2,0.1)
lines(itmb.new, itmb.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")


###############################################################################################
# itmc
#dev.new(width = 5, height = 5, units="in")
xlim=c(-35,-24)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="itm_c")

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(itmc[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}

res=cor.test(as.numeric(itmc), thin, method="pearson")
text(-32,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(-32,-28,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(itmc))
itmc.new=seq(-35,-24)
lines(itmc.new, itmc.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")


################################################################################
# enh
#dev.new(width = 5, height = 5, units="in")
xlim=c(1, 1.5)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab="enh")

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

   # plot
  points(enh[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(enh), thin, method="pearson")
text(1.1,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(1.1,-28,paste("p = ",signif(res$p.value, digits=4),sep=""), cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(enh))
enh.new=seq(1, 1.5,0.1)
lines(enh.new, enh.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")


########################################################################
# k
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,10)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(kappa))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(k[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(k), thin, method="pearson")
text(2,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(2,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(k))
k.new=seq(0,10)
lines(k.new, k.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")



############################################################################
# tau
#dev.new(width = 5, height = 5, units="in")
xlim=c(500,3000)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(tau))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(tau[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(tau), thin, method="pearson")
text(1000,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(1000,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(tau))
tau.new=seq(0,3000, 100)
lines(tau.new, tau.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")



##############################################################################
# lith
#dev.new(width = 5, height = 5, units="in")
xlim=c(50,200)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(lith))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(lith[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(lith), thin, method="pearson")
text(100,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(100,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(lith))
lith.new=seq(0,200,20)
lines(lith.new, lith.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")



###########################################################################
# kt
#dev.new(width = 5, height = 5, units="in")
xlim=c(0.0005,0.0015)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(kt))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(kt[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(kt), thin, method="pearson")
text(0.0008,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(0.0008,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(kt))
kt.new=seq(0.0005,0.0015,0.0001)
lines(kt.new, kt.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")


#####################################################################################
# fp
xlim=c(0.1,1.5)
ylim=c(-30,10)
plot(xlim, ylim, type="n", ylab="NGRIP thin rate (m/kyr)", xlab=expression(fp))

for(i in nr.best:1){
   #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  h = ncvar_get(nc,"H_ice")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  # calc thinning rat
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[length(time)]
  max.i=which(dzs.ngrip==max(dzs.ngrip))
  pd.i=length(time)
  min.i=pd.i
  thin.rate=(zs.ngrip[max.i]-zs.ngrip[min.i])/(time[max.i]-time[min.i])*1000  # m/kyr

  # plot
  points(fp[i],thin.rate, col="black",bg=mycol[i],pch=21, cex=1)
  print(i)
}
res=cor.test(as.numeric(fp), thin, method="pearson")
text(0.4,-25,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(0.4,-28,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin ~ as.numeric(fp))
fp.new=seq(0.1,1.5,0.1)
lines(fp.new, fp.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")




dev.off()
