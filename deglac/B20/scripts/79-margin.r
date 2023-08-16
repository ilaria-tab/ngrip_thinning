library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(viridis)

print("Start: 79_margin")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# NE transect-> 79N, Blaso
#xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698) # y coord
xs=c(200,744)
ys=c(-1354,-690)
fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)

# second transect
#xs=c(144, 848) # x coord for ZI GL position at the PD and max LGM extent (=cont. shelf break)
#ys=c(-1354,-874) # y coord
#fit=lm(ys~xs)
#x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
#y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
#dat=data.frame(x.line, y.line)

# GL from Larsen 2018
# Larsen curve averaged for all NE border
tab=read.table("/home/titan/gwgk/gwgk005h/work/ngrip/data/Larsen2018_margin_NEGIS_Fig2a-3a.dat",skip=0)
t.lar=tab$V1
gl.lar=tab$V2
#lines(t.lar, gl.lar, lty=2, lwd=1.5, col="#901388")

xlim=c(-20000,0)
ylim=c(-150,350)
mycol=rainbow(nr)
plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1.5)

rmse.gl=c()
for(run in 1:nr){
  #file = paste(out.fldr,"/itmb.",itmb[run],".itmc.",itmc[run],".btq.",q[run],".cbz0.",z0[run],".cfngs.",cfngs[run],".nffdlt.",neff[run],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",run-1,"/yelmo2D.nc", sep="")   

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.gl[run]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.gl[run]=NA
      next
    }
  }

  #mask=ncvar_get(nc, "mask_bed")
  mask=ncvar_get(nc, "f_ice")
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time=ncvar_get(nc,"time")
  #time=time-1950
  nc_close(nc)

  # where is the margin at the present?
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
    #if(mask[i,j,nbands]==4 | mask[i,j,nbands]==5){
     if(mask[i,j,nbands]==0){ 
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

      #if(mask[i,j,t]==4 | mask[i,j,nbands]==5){
      if(mask[i,j,t]==0){
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
  lines(time, d, col=mycol[run], lwd=2)

  # calc rmse between this gl distance and that from Larsen 2018
  t.i=which(time==-15000)
  t.f=length(time)
  mod=d[t.i:t.f]
  time=time[t.i:t.f]
  #mod=d

  obs=c()
  for(tt in 1:length(time)){
    ii=which(abs(time[tt]-t.lar)==min(abs(time[tt]-t.lar)))
    obs[tt]=gl.lar[ii]
  }

  diff = mod - obs
  sq.k = diff^2
  n = length(diff)

  rmse.gl[run]=sqrt(sum(as.vector(sq.k),na.rm=T)/n) 
  print(run)
}

lines(t.lar, gl.lar, lty=2, lwd=5, col="black")

file.out=paste(work.fldr,"/rmse-79-margin.txt",sep="")
write.table(rmse.gl, file.out, quote = F, sep = " ", row.names=F)


