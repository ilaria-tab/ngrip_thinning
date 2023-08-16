library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)
library(TeachingDemos)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test6")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test6")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")


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

###############################################################################################################


# One best sims
file = paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test6/itmb.-2.888.itmc.-22.064.btq.0.485.cbz0.-162.000.cfngs.0.179.nffdlt.0.049/yelmo2D.nc", sep="")
#file = paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test6/itmb.-2.676.itmc.-22.124.btq.0.975.cbz0.-225.000.cfngs.0.491.nffdlt.0.035/yelmo2D.nc", sep="")

nc=nc_open(file)
zs=ncvar_get(nc,"z_srf")
tsum=ncvar_get(nc,"Ta_sum")
tann=ncvar_get(nc,"Ta_ann")
pr=ncvar_get(nc,"Pr_ann")
dpr=ncvar_get(nc,"dPr_ann")
dtann=ncvar_get(nc,"dTa_ann")
smb=ncvar_get(nc,"smb")
us=ncvar_get(nc,"uxy_s")
h=ncvar_get(nc,"H_ice")
zb=ncvar_get(nc,"z_bed")
bm=ncvar_get(nc,"bmb_shlf")
visc=ncvar_get(nc,"visc_eff_int")
neff=ncvar_get(nc,"N_eff")
ty=ncvar_get(nc,"time")
ty=ty/1000
tsum=tsum-273.15
tann=tann-273.15
#nc_close(nc)
pd=nbands

zs.ngrip=zs[ngrip.i,ngrip.j,]
dzs.ngrip=zs.ngrip-zs.ngrip[pd]
tann.ngrip=tann[ngrip.i,ngrip.j,]
tsum.ngrip=tsum[ngrip.i,ngrip.j,]
pr.ngrip=pr[ngrip.i,ngrip.j,]
dpr.ngrip=dpr[ngrip.i,ngrip.j,]
smb.ngrip=smb[ngrip.i,ngrip.j,]
us.ngrip=us[ngrip.i,ngrip.j,]
h.ngrip=h[ngrip.i,ngrip.j,]
zb.ngrip=zb[ngrip.i,ngrip.j,]
visc.ngrip=visc[ngrip.i,ngrip.j,]
neff.ngrip=neff[ngrip.i,ngrip.j,]

en.i=egrip.i
en.j=egrip.j

zs.en=zs[en.i,en.j,]
dzs.en=zs.en-zs.en[pd]
tann.en=tann[en.i,en.j,]
tsum.en=tsum[en.i,en.j,]
pr.en=pr[en.i,en.j,]
smb.en=smb[en.i,en.j,]
us.en=us[en.i,en.j,]
h.en=h[en.i,en.j,]
zb.en=zb[en.i,en.j,]
visc.en=visc[en.i,en.j,]
neff.en=neff[en.i,en.j,]

dtann.neoff=dtann[169,307,]
tsum.neoff=tsum[169,307,]

#tsum.agas=tsum[agas.i, agas.j,]
#dtann.agas=dtann[agas.i, agas.j,]

# GL retreat
# third transect-> 79N, Blaso
xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
ys=c(-1322,-698) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)


  mask=ncvar_get(nc, "mask_bed")
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time=ncvar_get(nc,"time")

  # where is the GL at the present?
  for(pt in 1:dim(dat)[1]){
    px=round(dat[pt,1])
    py=round(dat[pt,2])
    # from coord to i,j
    i=which(abs(xc-px)==min(abs(xc-px)))
    j=which(abs(yc-py)==min(abs(yc-py)))
    if(length(i)>1 | length(j)>1){
      i=i[1]
      j=j[1]
    }
    if(mask[i,j,42]==4){
      gl.i=i
      gl.j=j
      break
    }else{
      # nothing
    }
  }
  pd.gl=c(gl.i,gl.j)
#  pd.gl=c(73,150)

  # calc distance from pd.gl
  d=c()
  for(t in 1:length(time)){
    for(pt in 1:dim(dat)[1]){
      px=round(dat[pt,1])
      py=round(dat[pt,2])

      # from coord to i,j
      i=which(abs(xc-px)==min(abs(xc-px)))
      j=which(abs(yc-py)==min(abs(yc-py)))

      if(length(i)>1 | length(j)>1){
        i=i[1]
        j=j[1]
      }

      if(mask[i,j,t]==4){
        gl.i=i
        gl.j=j
        break
      }else{
        # nothing
      }
    }
    d[t]=sqrt((xc[pd.gl[1]]-xc[gl.i])^2 + (yc[pd.gl[2]]-yc[gl.j])^2)
    if(gl.i<pd.gl[1] | gl.j<pd.gl[2]){d[t]=-d[t]}
  }

dtann.gl=dtann[gl.i,gl.j,]



plot.out=paste("/home/itabone/work/v1.23/lhs/B18/8km/test6/Figures/1D.png", sep="")
png(plot.out, width=5, height=10, units="in", res=100)
#dev.new(width=5, height=10, units="in", res=100)
#par(mfrow=c(5,1), mar=c(4,4,1.0,0.5), oma=c(0,0,0,0))
par(mfrow=c(6,1), mar=c(0,4,0,4), oma=c(4,3,1,4))

# S
at.x=seq(-15,0,by=2.5)
at.y=seq(-50,250, by=50)
xlim=c(-15,0)
ylim=c(-50,250)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
lines(t.ds.v09/1000, ds.ngrip.v09, col="grey80")
lines(t.ds.v09/1000, ds.ngrip.v09 + 40, col="grey80")
lines(t.ds.v09/1000, ds.ngrip.v09 - 40, col="grey80")
lines(t.ds.l13/1000, ds.ngrip1.l13, col="grey60",lty=2)
lines(t.ds.l13/1000, ds.ngrip2.l13, col="grey60", lty=2)
lines(t.ds.l13/1000, ds.ngrip.max.l13, col="grey60", lty=2)
lines(t.ds.l13/1000, ds.ngrip.min.l13, col="grey60", lty=2)
lines(ty, dzs.ngrip, col="#342E37", lwd=3)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("S (m)",side=2,line=4.5,cex=1)
#abline(v=-9.3, col="pink", lwd=2)

# U
at.y=seq(0,2, by=1)
ylim=c(0,2)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
lines(ty, log10(us.ngrip), lwd=3, col="#342E37")
lines(ty, log10(us.en), col="#6CC25B", lwd=3, lty=1)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("U log10(m/yr)",side=4,line=3.2,cex=1)
#abline(v=-9.3, col="pink", lwd=2)

# Neff
at.y=seq(7.4,7.5, by=0.05)
ylim=c(7.4,7.5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
lines(ty, log10(neff.en), lwd=3, col="#6CC25B")
#abline(v=-9.3, col="pink", lwd=2)
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Neff EGRIP log10(Pa)",side=2,line=5, cex=1)


# GL retreat
at.y=seq(-100,350, by=50)
ylim=c(-100,350)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#rect(-9.3,-70,-1,-20, col="mistyrose", border=NA)
lines(ty, d, col="#E93835", lwd=3)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("GL distance (km)",side=4,line=4, cex=1)
#abline(v=-9.3, col="pink", lwd=2)
# Data from Larsen 2018
abline(h=0, col="grey")
t.data=c(-7.90)
t.err=c(2.69)
for(i in 1:length(t.data)){
  arrows(x0=(t.data[i]-t.err[i]),y0=-50,x1=(t.data[i]+t.err[i]),y1=-50, code=3, angle=90, length=0.05, col="grey40",lwd=1)
}
points(t.data,-50, col="grey40")
# Data from Bennike Weidick 2001 (Blaso) but taken from Larsen2018 SM
t.min=c(-7441,-7252,-7156,-7151,-7150,-6970,-6930,-6393,-5928,-5941,-5566,-4820,-4790)
t.max=c(-7163,-6980,-6795,-6739,-6736,-6679,-6666,-6020,-5664,-5643,-5262,-4525,-4480)
t.min=t.min/1000
t.max=t.max/1000
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=-90,x1=t.max[i],y1=-90, code=3, angle=90, length=0.05, col="grey20", lwd=1)
    #arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
  }
# Larsen curve averaged for all NE border
tab=read.table("/mnt/lustre/scratch/itabone/ilaria/data/Greenland/Larsen2018_groundline_NEGIS.dat",skip=0)
t.lar=tab$V1/1000
gl.lar=tab$V2
lines(t.lar, gl.lar, lty=2, lwd=1.5, col="#901388")


# dTann 
file="/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18.nc"
nc=nc_open(file)
dtas.b18=ncvar_get(nc, "tas")
t.b18=ncvar_get(nc, "time")
nc_close(nc)
#dtas.b18.79=dtas.b18[156,306,,]
dtas.b18.79=dtas.b18[149,305,,]
#dtas.b18.79=dtas.b18[gl.i,gl.j,,]
dtann.b18.79=c()
for(t in 1:length(t.b18)){
  dtann.b18.79[t]=mean(dtas.b18.79[1:12,t])
}
# shift the anomaly wrt to 1750 CE (t[2181])
dtann.b18.79=dtann.b18.79-dtann.b18.79[2181]

at.y=seq(-10,7.5, by=2.5)
#xlim=c(-15,0)
ylim=c(-10,7.5)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
polygon(c(t.tann.l17/1000, rev(t.tann.l17/1000)), c(tann.agas.max.l17, rev(tann.agas.min.l17)), border = NA, col=alpha.col("grey50",40))
lines(t.tann.l17/1000,tann.agas.l17, lwd=2, col="grey50")
#lines(ty,dtann.neoff, lwd=3, col="#1098F7")
lines(t.b18/1000, dtann.b18.79, lwd=1.5, col="#1098F7")
axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Tann anomaly (°C)",side=2,line=4.5,cex=1)
#abline(v=-9.3, col="pink", lwd=2)



# Pr fraction
at.y=seq(0.2,1.2,by=0.2)
ylim=c(0.2,1.2)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
lines(ty, pr.ngrip/pr.ngrip[length(ty)], col="#342E37", lwd=3)
axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
mtext("Pr fraction",side=4,line=4, cex=1)
#abline(v=-9.3, col="pink", lwd=2)

axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1.3)
mtext("Time (kyr)",side=1,line=2.75,cex=1)


dev.off()

#####################################################

# plot map
file="/home/itabone/data/ice_data/Greenland/GRL-16KM/GRL-16KM_TOPO-M17.nc"
zs=raster(file, varname="z_srf")
zs[zs<1]=NA

plot.out=paste("/home/itabone/work/v1.23/lhs/B18/8km/test6/Figures/map-1D.png", sep="")
png(plot.out, width=4, height=5, units="in", res=100)
#dev.new(width=4, height=5, units="in")
par(mar=c(1,1,1,1))
plot(zs, col=c("#dee2e6"), box=F,axes=F,legend=F)
contour(zs, add=T, drawlabels=F)
points(xc[ngrip.i],yc[ngrip.j],col="#342E37",pch=20, cex=2)
points(xc[egrip.i],yc[egrip.j],col="#6CC25B", pch=20, cex=2)
points(xc[gl.i],yc[gl.j],col="#E93835", pch=20, cex=2)
points(xc[156],yc[306],col="#1098F7", pch=20, cex=2)
dev.off()


###################################################
# caveats

# Tsum
#at.y=seq(-10,10, by=5)
#ylim=c(-10,10)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#lines(ty,tsum.neoff, lwd=3, col="green")
#lines(t.tsum.l17/1000,tsum.agas.l17, lwd=2, col="magenta")
#axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("Tsum (°C)",side=4,line=3.2,cex=1)
#abline(v=-9.3, col="pink", lwd=2)

## Tsum from B18
#file="/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18.nc"
#nc=nc_open(file)
#dtas.b18=ncvar_get(nc, "tas")
#t.b18=ncvar_get(nc, "time")
#nc_close(nc)
#
#dtas.b18.79=dtas.b18[147,300,,]
#dtsum.b18.79=c()
#for(t in 1:length(t.b18)){
#    dtsum.b18.79[t]=mean(dtas.b18.79[6:8,t])
#}


## Pr
#at.y=seq(0.02,0.22, by=0.05)
#ylim=c(0.02,0.22)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#lines(t.k14/1000, acc.ngrip.k14, col="grey40")
#lines(t.ngrip.b20/1000, acc.ngrip.b20, col="grey40", lty=2)
#lines(t.g14/1000, acc.ngrip.g14, col="grey40",lty=3)
#lines(ty, pr.ngrip, col="black", lwd=3)
#legend("topleft", c("k14","b20","g14"), lty=c(1,2,3),bty="n", cex=1.5)
#axis(4,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("Accum (m/yr)",side=4,line=4, cex=1)
#abline(v=-9.3, col="pink", lwd=2)

# Visc
#at.y=seq(1e11,3e12, by=6e11)
#ylim=c(1e11,3e12)
#plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)
#lines(ty, visc.ngrip, lwd=3)
#abline(v=-9.3, col="orange", lwd=2)
#axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1.3)
#mtext("Visc (Pa m a)",side=2,line=5, cex=1)


