library(fields)
#library(raster)
library(ncdf4)

# ice core locations
filename = paste("~/ilaria/scripts/ice-core-16km.txt")
txt=read.table(filename,header=F, sep="", skip=3)
i=txt$V2
j=txt$V3
ngrip.i=i[1]
ngrip.j=j[1]
grip.i=i[2]
grip.j=j[2]
dye3.i=i[3]
dye3.j=j[3]
campcen.i=i[4]
campcen.j=j[4]

ngrip.zs=2920

q=c("0.0","0.3","0.6","1.0")
for(i in q){
   file=paste("~/ilaria/outputs/btq.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.q.",i,sep=""),dzs)
   nc_close(nc)
}

cf=c("0.25","0.5","0.75","1.0")
for(i in cf){
   file=paste("~/ilaria/outputs/cfstrm.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.cf.",i,sep=""),dzs)
   nc_close(nc)
}

enhshr=c("1","3","5")
for(i in enhshr){
   file=paste("~/ilaria/outputs/enhshr.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.enhshr.",i,sep=""),dzs)
   nc_close(nc)
}

file="~/ilaria/outputs/ghf-M18/yelmo2D.nc"
nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
time=ncvar_get(nc,"time")
time=time/1000
zs.core=zs[ngrip.i,ngrip.j,]
dzs.ghf=zs.core-zs.core[pd]
nc_close(nc)

itmb=c("-5","-3","-1","0")
for(i in itmb){
   file=paste("~/ilaria/outputs/itmb.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.itmb.",-as.numeric(i),sep=""),dzs)
   nc_close(nc)
}

itmc=c("-15","-20","-30")
for(i in itmc){
   file=paste("~/ilaria/outputs/itmc.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.itmc.",-as.numeric(i),sep=""),dzs)
   nc_close(nc)
}

kppgrz=c(10,50,100,200,500,1000)
for(i in kppgrz){
   file=paste("~/ilaria/outputs/kppgrz.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.kppgrz.",i,sep=""),dzs)
   nc_close(nc)
}

nffdlt=c(0.02,0.03)
for(i in nffdlt){
   file=paste("~/ilaria/outputs/nffdlt.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.nffdlt.",i,sep=""),dzs)
   nc_close(nc)
}

tau=c(1000,2000,5000)
for(i in tau){
   file=paste("~/ilaria/outputs/t.",i,"/yelmo2D.nc",sep="")
   nc=nc_open(file)
   zs=ncvar_get(nc,"z_srf")
   time=ncvar_get(nc,"time")
   time=time/1000
   zs.core=zs[ngrip.i,ngrip.j,]
   pd=length(time)
   dzs=zs.core-zs.core[pd]
   assign(paste("dzs.tau.",i,sep=""),dzs)
   nc_close(nc)
}


plot(time, dzs.q.0.0, xlim=c(-10, 0), ylim=c(-10,85), type="l",lty=1,col="grey")
lines(time, dzs.q.0.3, col="grey", lty=1)
lines(time, dzs.q.0.6, col="grey", lty=2)
lines(time, dzs.q.1.0, col="grey", lty=3)
lines(time, dzs.ghf, col="red", lty=1)
lines(time, dzs.cf.0.25, col="orange", lty=1)
lines(time, dzs.cf.0.5, col="orange", lty=2)
lines(time, dzs.cf.0.75, col="orange", lty=3)
lines(time, dzs.cf.1.0, col="orange", lty=4)
lines(time, dzs.enhshr.1, col="green", lty=1)
lines(time, dzs.enhshr.3, col="green", lty=2)
lines(time, dzs.enhshr.5, col="green", lty=3)
lines(time, dzs.itmb.1, col="yellow", lty=1)
lines(time, dzs.itmb.3, col="yellow", lty=2)
lines(time, dzs.itmb.5, col="yellow", lty=3)
lines(time, dzs.itmc.15, col="lightblue", lty=1)
lines(time, dzs.itmc.20, col="lightblue", lty=2)
lines(time, dzs.itmc.30, col="lightblue", lty=3)
lines(time, dzs.kppgrz.10, col="magenta", lty=1)
lines(time, dzs.kppgrz.100, col="magenta", lty=2)
lines(time, dzs.kppgrz.500, col="magenta", lty=3)
lines(time, dzs.kppgrz.1000, col="magenta", lty=4)
lines(time, dzs.nffdlt.0.02, col="brown", lty=1)
lines(time, dzs.nffdlt.0.03, col="brown", lty=2)
lines(time, dzs.tau.1000, col="blue", lty=1)
lines(time, dzs.tau.2000, col="blue", lty=1)


