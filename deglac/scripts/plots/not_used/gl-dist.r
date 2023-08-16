# Plots

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


##########################################################################################################################
# define transect 
xs=c(240,753) # x coord for ZI GL position at the PD and max LGM extent (=cont. shelf break)
ys=c(-1546,-698) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)

#plot(zs)
#segments(dat[1,1],dat[1,2],dat[dim(dat)[1],1], dat[dim(dat)[1],2], lty=1, lwd=3, col="deeppink")

# second transect
xs=c(144, 848) # x coord for ZI GL position at the PD and max LGM extent (=cont. shelf break)
ys=c(-1354,-874) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)

# third transect-> 79N, Blaso
xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
ys=c(-1322,-698) # y coord
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)


# plot distance vs time
plot.out=paste(work.fldr,"/Figures/gldist-vs-time.png", sep="")
png(plot.out, width = 7, height = 5, units = "in", res=100)
#dev.new(width = 7, height = 5, units="in")
xlim=c(-20000,0)
ylim=c(-200,350)
plot(xlim, ylim, type="n", xlab="Time (yr)", ylab="GL distance (km)")

for(run in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[run],".itmc.",itmc[run],".btq.",q[run],".cbz0.",z0[run],".cfngs.",cfngs[run],".nffdlt.",neff[run],"/yelmo2D.nc", sep="")

  nc=nc_open(file)
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
    if(mask[i,j,length(time)]==4){
      gl.i=i
      gl.j=j
      break
    }else{
      # nothing
    }
  }
  pd.gl=c(gl.i,gl.j)
#  pd.gl=c(145,300)

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
  # plot
  lines(time, d, col=mycol[run], lwd=1)
}
abline(h=0, col="black")

# Data from Larsen 2018
#abline(h=0, col="grey")
t.data=c(-7900)
t.err=c(2690)
for(i in 1:length(t.data)){
  arrows(x0=(t.data[i]-t.err[i]),y0=-50,x1=(t.data[i]+t.err[i]),y1=-50, code=3, angle=90, length=0.05, col="grey40",lwd=1)
}
points(t.data,-50, col="grey40")

# Data from Bennike Weidick 2001 (Blaso) but taken from Larsen2018 SM
t.min=c(-7441,-7252,-7156,-7151,-7150,-6970,-6930,-6393,-5928,-5941,-5566,-4820,-4790)
t.max=c(-7163,-6980,-6795,-6739,-6736,-6679,-6666,-6020,-5664,-5643,-5262,-4525,-4480)
t.min=t.min
t.max=t.max
for(i in 1:length(t.max)){
    arrows(x0=t.min[i],y0=-150,x1=t.max[i],y1=-150, code=3, angle=90, length=0.05, col="grey20", lwd=1)
    #arrows(x0=t.min[i]+(t.max[i]-t.min[i])/2,y0=rsl[i]-rsl.err[i],x1=t.max[i]-(t.max[i]-t.min[i])/2,y1=rsl[i]+rsl.err[i], code=3, angle=90,length=0.05, col="red", lwd=1)
  }

# Larsen curve averaged for all NE border
tab=read.table("/mnt/lustre/scratch/itabone/ilaria/data/Greenland/Larsen2018_groundline_NEGIS.dat",skip=0)
t.lar=tab$V1
gl.lar=tab$V2
lines(t.lar, gl.lar, lty=2, lwd=4, col="black")

dev.off()

