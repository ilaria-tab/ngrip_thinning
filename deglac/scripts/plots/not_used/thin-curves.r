library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

# col palette
library(viridis)
#mycol=tim.colors(nr)
col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
mycol=col[as.numeric(cut(b.best[,9],breaks = nr))]

alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}


# GRIP
plot.out=paste(work.fldr,"/Figures/B20-elev-grip.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in")
par(mfrow=c(2,1), mar=c(4,4,1.0,0.5), oma=c(0,0,0,0)) 
# dS
plot(t.ds.v09/1000, ds.grip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="GRIP", ylim=c(-50,250), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.grip.v09 - 40, rev(ds.grip.v09 +40)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.grip.min.l13, rev(ds.grip.max.l13)), border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, ds.grip1.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, ds.grip2.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.grip=zs[grip.i,grip.j,]
  dzs.grip=zs.grip-zs.grip[pd]  

  #if(b.ord[i,9] > 100){
    lines(time,  dzs.grip, col=alpha.col(mycol[i],100), lty=1, lwd=2)
  #}else{
  #  next
  #}
}

# S
plot(t.ds.v09/1000, s.grip.v09, type="l", xlab="Time (kyr)", ylab="S (m)",main="GRIP", ylim=c(3100,3500), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.grip.v09 - 40, rev(s.grip.v09 +40)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.grip.min.l13, rev(s.grip.max.l13)),  border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, s.grip1.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, s.grip2.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.grip=zs[grip.i,grip.j,]
  dzs.grip=zs.grip-zs.grip[pd]

  #if(b.ord[i,9] > 100){
    lines(time,  zs.grip, col=alpha.col(mycol[i],100), lty=1, lwd=2)
  #}else{
  #  next
  #}
}
dev.off()



# NGRIP
plot.out=paste(work.fldr,"/Figures/B20-elev-ngrip.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in")
par(mfrow=c(2,1), mar=c(4,4,1.0,0.5), oma=c(0,0,0,0))
# dS
plot(t.ds.v09/1000, ds.ngrip.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="NGRIP", ylim=c(-50,250), xlim=c(-15,0),col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.ngrip.v09 - 40, rev(ds.ngrip.v09 +40)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.ngrip.min.l13, rev(ds.ngrip.max.l13)), border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, ds.ngrip1.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, ds.ngrip2.l13, col="grey70", lty=1, lwd=2)
#nr=15
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")

  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
 
  #if(b.ord[i,9] > 100){
    lines(time,  dzs.ngrip, col=alpha.col(mycol[i],100), lty=1, lwd=2)
  #}
}

# S
plot(t.ds.v09/1000, s.ngrip.v09, type="l", xlab="Time (kyr)", ylab="S (m)",main="NGRIP", ylim=c(2800,3200),xlim=c(-15,0), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.ngrip.v09 - 40, rev(s.ngrip.v09 +40)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.ngrip.min.l13, rev(s.ngrip.max.l13)),  border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, s.ngrip1.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, s.ngrip2.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  dzs.ngrip=zs.ngrip-zs.ngrip[pd]
  
  #if(b.ord[i,9] > 100){
    lines(time,  zs.ngrip, col=alpha.col(mycol[i],100), lty=1, lwd=2)
  #}else{
  # next
  #}
}
dev.off()


# DYE3
plot.out=paste(work.fldr,"/Figures/B20-elev-dye3.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in")
par(mfrow=c(2,1), mar=c(4,4,1.0,0.5), oma=c(0,0,0,0))
# dS
plot(t.ds.v09/1000, ds.dye3.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="DYE3", ylim=c(-150,500), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.dye3.v09 - 65, rev(ds.dye3.v09 +65)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(ds.dye3.min.l13, rev(ds.dye3.max.l13)), border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, ds.dye31.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, ds.dye32.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".kppgrz.",kgrz[i],".btq.",q[i],".cbz0.",z0[i],".cfn.",cfn[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.dye3=zs[dye3.i,dye3.j,]
  dzs.dye3=zs.dye3-zs.dye3[pd]
  lines(time,  dzs.dye3, col=alpha.col(mycol[1],100), lty=1, lwd=2)
}
# S
plot(t.ds.v09/1000, s.dye3.v09, type="l", xlab="Time (kyr)", ylab="S (m)",main="DYE3", ylim=c(2300,3000), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.dye3.v09 - 65, rev(s.dye3.v09 +65)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l13/1000, rev(t.ds.l13/1000)), c(s.dye3.min.l13, rev(s.dye3.max.l13)),  border = NA, col=alpha.col("grey80",50))
lines(t.ds.l13/1000, s.dye31.l13, col="grey70",lty=1, lwd=2)
lines(t.ds.l13/1000, s.dye32.l13, col="grey70", lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".kppgrz.",kgrz[i],".btq.",q[i],".cbz0.",z0[i],".cfn.",cfn[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.dye3=zs[dye3.i,dye3.j,]
  dzs.dye3=zs.dye3-zs.dye3[pd]
  lines(time,  zs.dye3, col=alpha.col(mycol[1],100), lty=1, lwd=2)
}
dev.off()



# CAMP CENTURY
plot.out=paste(work.fldr,"/Figures/B20-elev-campcen.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in")
par(mfrow=c(2,1), mar=c(4,4,1.0,0.5), oma=c(0,0,0,0))
# dS
plot(t.ds.v09/1000, ds.campcen.v09, type="l", xlab="Time (kyr)", ylab="dS (m)",main="CAMPCEN", ylim=c(-200,1100), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(ds.campcen.v09 - 65, rev(ds.campcen.v09 +65)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l17/1000, rev(t.ds.l17/1000)), c(ds.campcen.min.l17, rev(ds.campcen.max.l17)), border = NA, col=alpha.col("grey80",50))
lines(t.ds.l17/1000, ds.campcen.l17, col="grey70",lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".kppgrz.",kgrz[i],".btq.",q[i],".cbz0.",z0[i],".cfn.",cfn[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.campcen=zs[campcen.i,campcen.j,]
  dzs.campcen=zs.campcen-zs.campcen[pd]
  lines(time,  dzs.campcen, col=alpha.col(mycol[1],100), lty=1, lwd=2)
}
# S
plot(t.ds.v09/1000, s.campcen.v09, type="l", xlab="Time (kyr)", ylab="S (m)",main="CAMPCEN", ylim=c(1700,3000), col="#9FC8D2", lwd=2)
polygon(c(t.ds.v09/1000, rev(t.ds.v09/1000)), c(s.campcen.v09 - 65, rev(s.campcen.v09 +65)), border = NA, col=alpha.col("#9FC8D2", 50))
polygon(c(t.ds.l17/1000, rev(t.ds.l17/1000)), c(s.campcen.min.l17, rev(s.campcen.max.l17)),  border = NA, col=alpha.col("grey80",50))
lines(t.ds.l17/1000, s.campcen.l17, col="grey70",lty=1, lwd=2)
for(i in 1:nr){
  file = paste(out.fldr,"/itmb.",itmb[i],".kppgrz.",kgrz[i],".btq.",q[i],".cbz0.",z0[i],".cfn.",cfn[i],"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  time=time/1000
  pd=length(time)
  nc_close(nc)

  zs.campcen=zs[campcen.i,campcen.j,]
  dzs.campcen=zs.campcen-zs.campcen[pd]
  lines(time,  zs.campcen, col=alpha.col(mycol[1],100), lty=1, lwd=2)
}
dev.off()



























