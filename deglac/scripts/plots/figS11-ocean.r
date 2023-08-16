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

#######################################################################################################
### ocean
i=183
j=311

file=paste(out.fldr,"/",b.best[1,1]-1,"/yelmo2D.nc", sep="")
nc = nc_open(file)
xc=ncvar_get(nc, "xc")
yc= ncvar_get(nc, "yc")
time = ncvar_get(nc,"time")
bmb = ncvar_get(nc, "bmb_shlf")
dta=ncvar_get(nc,"dTa_ann")
dto=ncvar_get(nc,"dT_shlf")
nc_close(nc)
	
file=paste(out.fldr,"/tests/trace21/",b.best[1,1]-1,"/yelmo2D.nc", sep="")
nc = nc_open(file)
time.trace = ncvar_get(nc,"time")
dto.trace=ncvar_get(nc,"dT_shlf")
bmb.trace=ncvar_get(nc, "bmb_shlf")
nc_close(nc)

bmb=bmb[i,j,]
dta=dta[i,j,]
dto=dto[i,j,]
dto.trace=dto.trace[i,j,]
bmb.trace=bmb.trace[i,j,]
time=time/1000
time.trace=time.trace/1000

plot.out=paste(work.fldr,"/Figures/figS11-ocean.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in", res=100)
par(mar=c(4,4,0,4), oma=c(0,3,1,2), mfrow=c(3,1))

plot(time, dta,  xlim=c(-15,0), type="l",lwd=2, xlab="", ylab="Annual Tair anomaly (°C)")
grid()
plot(time, dto,  xlim=c(-15,0), type="l", lwd=2, xlab="",ylab="Annual Tocn anomaly (°C)")
lines(time.trace, dto.trace, col="blue", lwd=2)
grid()
plot(time, bmb, xlim=c(-15,0), type="l", lwd=2,xlab="Time (kyr)",ylab="Submarine melt (m/yr)")
lines(time.trace, bmb.trace, col="blue", lwd=2)
grid()


dev.off()


