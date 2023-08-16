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

# col palette runs based on skill score
#col=colorRampPalette(c("blue","lightblue","green","yellow","orange","red","darkred"))(nr)
col=rev(plasma(nr.best))
mycol=col[as.numeric(cut(b.best[,17],breaks = nr.best))]

# col alpha
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_alpha_col.r")


###############################################################################################################
file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
#time=time-1950
#nc_close(nc)
  #nc = nc_open(file)
ub = ncvar_get(nc, "uxy_b")
h =  ncvar_get(nc, "H_ice")
hw =  ncvar_get(nc, "H_w")
neff =  ncvar_get(nc, "N_eff")
bmb =  ncvar_get(nc, "bmb")
beta =  ncvar_get(nc, "beta")

i=139
j=271

ub=ub[i,j,]
h= h[i,j,]
hw =hw[i,j,]
neff =neff[i,j,]
bmb=bmb[i,j,]
beta=beta[i,j,]
nc_close(nc)
time=time/1000

plot.out=paste(work.fldr,"/Figures/figS12-dynamics.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in", res=100)
#layout(matrix(c(1,2,3),nrow=3,byrow=T), heights=c(1,1,1))
#layout.show(n=6)
par(mar=c(4,4,0,4), oma=c(0,3,1,2), mfrow=c(6,1))

plot(time, h,  xlim=c(-15,0), type="l", xlab="", ylab="H ice (m)")
grid()
plot(time, neff,  xlim=c(-15,0), type="l", xlab="",ylab="Neff (Pa)")
grid()
plot(time, beta, xlim=c(-15,0), type="l", xlab="",ylab="Beta (Pa yr m-1)")
grid()
plot(time, hw,  xlim=c(-15,0), type="l", xlab="",ylab="H water (m)")
grid()
plot(time, bmb,  xlim=c(-15,0), type="l",xlab="", ylab="Bmb (m/yr)")
grid()
plot(time, ub,  xlim=c(-15,0), type="l", xlab="Time (kyr)", ylab="Ub (m/yr)")
grid()
  

dev.off()

