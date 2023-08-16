library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/test45_deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/test45_deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
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

###############################################################################################################

i=1

file=paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc",sep="")
nc = nc_open(file)
sl=ncvar_get(nc,"z_sl")
zb=ncvar_get(nc,"z_bed")
time=ncvar_get(nc,"time")
qx=ncvar_get(nc,"qq_acx") # 0.5(H_i+H_i+1)*dx*ux (m3/a)
qy=ncvar_get(nc,"qq_acy")
smb=ncvar_get(nc,"smb")
bmb=ncvar_get(nc,"bmb")
hice=ncvar_get(nc,"H_ice")
taub=ncvar_get(nc,"taub")
nc_close(nc)

taub.egrip=taub[egrip.i,egrip.j,]
taub.ngrip=taub[ngrip.i,ngrip.j,]
hice.ngrip=hice[ngrip.i,ngrip.j,]
smb.ngrip=smb[ngrip.i,ngrip.j,]   #m/a
qx.ngrip=qx[ngrip.i,ngrip.j,]/8000/8000 #m/a
qy.ngrip=qy[ngrip.i,ngrip.j,]/8000/8000 # m/a
q.ngrip=sqrt((qx.ngrip)^2+(qy.ngrip)^2) # m/a
#q.ngrip=q.ngrip/8000/8000 # m/a
bmb.ngrip=bmb[ngrip.i,ngrip.j,] 

dhdt=c()
for(t in 1:length(time)){
   dhdt[t]=(hice.ngrip[t]-hice.ngrip[t-1])/(time[t]-time[t-1])
}

dev.new(width=5, height=8, units="in", res=100)
layout(matrix(c(1,2,3,4),nrow=4,byrow=T), heights=c(1,1,1,1))

par(mar=c(0,4,0,4), oma=c(4,3,1,2))
#plot(time, taub.egrip, type="l", xlim=c(-15000, -6000))
#abline(v=-12200, col="green")
plot(time, hice.ngrip, type="l", xlim=c(-15000, -6000))
abline(v=-11000, col="red")
abline(v=-9800, col="magenta")
abline(v=-10600, col="orange")
#plot(time, qx.ngrip+qy.ngrip, type="l", xlim=c(-15000, -6000))
#abline(v=-11200, col="yellow")
plot(time, q.ngrip, type="l", xlim=c(-15000, -6000))
abline(v=-11000, col="red")
abline(v=-9800, col="magenta")
abline(v=-10600, col="orange")
plot(time, smb.ngrip+bmb.ngrip, type="l", xlim=c(-15000, -6000))
abline(v=-11000, col="red")
abline(v=-9800, col="magenta")
abline(v=-10600, col="orange")
plot(time, dhdt, type="l", xlim=c(-15000, -6000))
abline(v=-11000, col="red")
abline(v=-9800, col="magenta")
abline(v=-10600, col="orange")
abline(h=0, lty=3)



#dev.new()
#plot(time, smb.ngrip-q.ngrip, type="l")
#abline(h=0, lty=3)
#
## rate of change
#for(t in 1:(length(time)-1)){
#  dhdt[t]=(hice.ngrip[t+1]-hice.ngrip[t])/(time[t+1]-time[t])
#}
#dhdt[90]=dhdt[89]
#
#dev.new()
#plot(time,dhdt, type="l")
#abline(h=0)
#
#grid()
