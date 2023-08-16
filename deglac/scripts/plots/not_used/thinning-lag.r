# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)
library(TeachingDemos)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/B18-TAS-B20-highP")
work.fldr=paste("/home/itabone/work/ngrip/v1.23/lhs/B18/8km/B18-TAS-B20-highP")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/ngrip/data/load-ice-core-data-8km.r")


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


##########################################################################################################################
#source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))
#t.zmax=t.taumax=t.neffmax=t.dqbmax=t.dticemax=c()
#t.zmax.best=t.taumax.best=t.neffmax.best=t.dqbmax.best=t.dticemax.best=c()
#for(i in 1:nr){
#  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
#  nc=nc_open(file)
#  zs=ncvar_get(nc,"z_srf")
#  beta=ncvar_get(nc,"beta")
#  taub=ncvar_get(nc,"taub")
#  n_eff=ncvar_get(nc,"N_eff")
#  qb=ncvar_get(nc,"Q_b")
#  tice=ncvar_get(nc,"T_ice")
#  tice=tice[,,1,] # ice base
#  basins=ncvar_get(nc,"basins")
#  ty=ncvar_get(nc,"time")
#  ty=ty/1000
#
#  zs.ngrip=zs[ngrip.i,ngrip.j,]
#
#  taub.negis=neff.negis=qb.negis=tice.negis=c()
#  for(t in 1:length(ty)){
#    taub.negis[t]=mean(taub[,,t][basins==9])
#    neff.negis[t]=mean(n_eff[,,t][basins==9])
#    qb.negis[t]=mean(qb[,,t][basins==9])
#    tice.negis[t]=mean(tice[,,t][basins==9])
#  }
#
#  i.zmax=which(zs.ngrip==max(zs.ngrip))
#  i.taumax=which(taub.negis==max(taub.negis))
#  i.neffmax=which(neff.negis==max(neff.negis))
#
#  dqb.negis=c()
#  dtice.negis=c()
#  for(l in 1:(length(ty)) ){
#     dqb.negis[l]=(qb.negis[1+l]-qb.negis[l])/(ty[l+1]-ty[l])
#     dtice.negis[l]=(tice.negis[1+l]-tice.negis[l])/(ty[l+1]-ty[l])
#  }
#  i.dqbmax=which(dqb.negis==max(dqb.negis,na.rm=T))
#  dtice.negis[1:4]=dtice.negis[5]
#  i.dticemax=which(dtice.negis==max(dtice.negis, na.rm=T))
# 
#  t.zmax[i]=ty[i.zmax]
#  t.taumax[i]=ty[i.taumax]
#  t.neffmax[i]=ty[i.neffmax] 
#  t.dqbmax[i]=ty[i.dqbmax]
#  t.dticemax[i]=ty[i.dticemax]
#}
#
## write results to a file
#file.out=paste(work.fldr,"/scripts/plots/tzmax.txt",sep="")
#write.table(t.zmax, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/ttaumax.txt",sep="")
#write.table(t.taumax, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tneffmax.txt",sep="")
#write.table(t.neffmax, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tdqbmax.txt",sep="")
#write.table(t.dqbmax, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tdticemax.txt",sep="")
#write.table(t.dticemax, file.out, quote = F, sep = " ", row.names=F)
#
#source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
#nr.best=nrow(b.best)
#for(i in 1:nr.best){
#  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
#  nc=nc_open(file)
#  zs=ncvar_get(nc,"z_srf")
#  beta=ncvar_get(nc,"beta")
#  taub=ncvar_get(nc,"taub")
#  n_eff=ncvar_get(nc,"N_eff")
#  qb=ncvar_get(nc,"Q_b")
#  tice=ncvar_get(nc,"T_ice")
#  tice=tice[,,1,] # ice base
#  basins=ncvar_get(nc,"basins")
#  ty=ncvar_get(nc,"time")
#  ty=ty/1000
#
#  zs.ngrip=zs[ngrip.i,ngrip.j,]
#
#  taub.negis=neff.negis=qb.negis=tice.negis=c()
#  for(t in 1:length(ty)){
#    taub.negis[t]=mean(taub[,,t][basins==9])
#    neff.negis[t]=mean(n_eff[,,t][basins==9])
#    qb.negis[t]=mean(qb[,,t][basins==9])
#    tice.negis[t]=mean(tice[,,t][basins==9])
#  }
#
#  i.zmax=which(zs.ngrip==max(zs.ngrip))
#  i.taumax=which(taub.negis==max(taub.negis))
#  i.neffmax=which(neff.negis==max(neff.negis))
#
#  dqb.negis=c()
#  dtice.negis=c()
#  for(l in 1:(length(ty)) ){
#     dqb.negis[l]=(qb.negis[1+l]-qb.negis[l])/(ty[l+1]-ty[l])
#     dtice.negis[l]=(tice.negis[1+l]-tice.negis[l])/(ty[l+1]-ty[l])
#  }
#  i.dqbmax=which(dqb.negis==max(dqb.negis,na.rm=T))
#  dtice.negis[1:4]=dtice.negis[5]
#  i.dticemax=which(dtice.negis==max(dtice.negis, na.rm=T))
#
#  t.zmax.best[i]=ty[i.zmax]
#  zmax.best[i]=zmax(i.zmax)
#  t.taumax.best[i]=ty[i.taumax]
#  t.neffmax.best[i]=ty[i.neffmax]
#  t.dqbmax.best[i]=ty[i.dqbmax]
#  t.dticemax.best[i]=ty[i.dticemax]
#}
#
## write results to a file
#file.out=paste(work.fldr,"/scripts/plots/tzmax-best.txt",sep="")
#write.table(t.zmax.best, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/ttaumax-best.txt",sep="")
#write.table(t.taumax.best, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tneffmax-best.txt",sep="")
#write.table(t.neffmax.best, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tdqbmax-best.txt",sep="")
#write.table(t.dqbmax.best, file.out, quote = F, sep = " ", row.names=F)
#file.out=paste(work.fldr,"/scripts/plots/tdticemax-best.txt",sep="")
#write.table(t.dticemax.best, file.out, quote = F, sep = " ", row.names=F)
#
#
#################################################################################################
t.zmax.best=read.table(paste(work.fldr,"/scripts/plots/tzmax-best.txt",sep=""), skip=1)
t.zmax=read.table(paste(work.fldr,"/scripts/plots/tzmax.txt",sep=""), skip=1)
t.taumax.best=read.table(paste(work.fldr,"/scripts/plots/ttaumax-best.txt",sep=""), skip=1)
t.taumax=read.table(paste(work.fldr,"/scripts/plots/ttaumax.txt",sep=""), skip=1)
t.dqbmax.best=read.table(paste(work.fldr,"/scripts/plots/tdqbmax-best.txt",sep=""), skip=1)
t.dqbmax=read.table(paste(work.fldr,"/scripts/plots/tdqbmax.txt",sep=""), skip=1)
t.dticemax.best=read.table(paste(work.fldr,"/scripts/plots/tdticemax-best.txt",sep=""), skip=1)
t.neffmax.best=read.table(paste(work.fldr,"/scripts/plots/tneffmax-best.txt",sep=""), skip=1)

t.zmax.best=unlist(t.zmax.best, use.names=FALSE)
t.zmax=unlist(t.zmax, use.names=FALSE)
t.taumax.best=unlist(t.taumax.best, use.names=FALSE)
t.taumax=unlist(t.taumax, use.names=FALSE)
t.dqbmax.best=unlist(t.dqbmax.best, use.names=FALSE)
t.dticemax.best=unlist(t.dticemax.best, use.names=FALSE)
t.neffmax.best=unlist(t.neffmax.best, use.names=FALSE)

lagtau.best=t.zmax.best-t.taumax.best
lagdtice.best=t.zmax.best-t.dticemax.best
lagdqb.best=t.zmax.best-t.dqbmax.best
lagneff.best=t.zmax.best-t.neffmax.best
#lag=t.zmax-t.taumax
#lag[lag>20]=NA
lagtau.best[lagtau.best>20]=NA

hist(lagtau.best, breaks=100, prob=T, xlim=c(-5,5),col="blue", xlab="Lag between onset of NEGIS formation and NGRIP thinning (kyr)")
lines(density(lagtau.best, na.rm=T), col="blue", lwd=2)



##p1=hist(lagtau, breaks=100)
#p2=hist(lagtau.best, breaks=500)
#p3=hist(lagdtice.best, breaks=100)
#
dens=density(lagtau.best, bw = 0.2, na.rm=T)
plot(dens, frame = FALSE, col = "steelblue", xlim=c(-2,3), main="",xlab="Lag between onset of NEGIS formation and NGRIP thinning (kyr)") 
polygon(dens, col = "steelblue")

dens=density(lagdtice.best,bw=0.2, na.rm=T)
lines(dens, col = "green")
polygon(dens, col = alpha.col("green",50))

dens=density(lagdqb.best, bw=0.2, na.rm=T)
lines(dens, col = "pink")
polygon(dens, col = alpha.col("pink",50))

abline(v=0, lty=2, lwd=2, col="grey")

legend("topright", c("Taub","dTice(z0)","dQb"), col=c("steelblue",alpha.col("green",50),alpha.col("pink",50)), lwd=4)
#dens=density(lagneff.best, na.rm=T)
#lines(dens, col = "yellow")
#polygon(dens, col = alpha.col("yellow",50))

##scale=max(p2$counts)/max(p2$density)
#
#plot(p2, xlim=c(-2,2),col="blue", ylim=c(0,2), xlab="Lag between onset of NEGIS formation and NGRIP thinning (kyr)")
#plot(p3, xlim=c(-2,2),col=alpha.col("green",40), add=T)
#abline(v=0, lty=2)
##lines(list(x = dens$x, y = scale * dens$y), col = "red", lwd = 2)
######################################################################################





