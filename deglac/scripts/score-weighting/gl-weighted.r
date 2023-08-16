# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
#library(plot3D)
library(viridis)
#library(TeachingDemos)

print("grounding line")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)
ss.norm=b.best[,17]

# load surface elevation and climatological data at ice core sites
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")

# GL retreat
# third transect-> 79N, Blaso
xs=c(288, 640) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
ys=c(-1322,-698) # y coord
# new segment
#xs=c(200,744)
#ys=c(-1354,-690)
#xs=c(304,640)
#ys=c(-1322,-698)
#xs=c(288, 704) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698)  # good!
#xs=c(272, 704) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698)
#xs=c(264, 704) # x coord for GL extended to the interior  and max LGM extent (=cont. shelf break)
#ys=c(-1322,-698)  # good!

fit=lm(ys~xs)
x.line=seq(xs[1],xs[2])  # take into account any retreat (advance) from PD (LGM) GL position
y.line=fit$coefficients[[1]] + x.line*fit$coefficients[[2]]
dat=data.frame(x.line, y.line)

# score-weighted mean
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="") 

  # does the file exist?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    nom[i,]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      nom[i,]=NA
      next
    }
  }


#nc = nc_open(file)
  mask=ncvar_get(nc, "f_ice")
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time=ncvar_get(nc,"time")
  nc_close(nc)
 
  # where is the GL at the present?
  for(pt in 1:dim(dat)[1]){
    px=round(dat[pt,1])
    py=round(dat[pt,2])
    # from coord to i,j
    ii=which(abs(xc-px)==min(abs(xc-px)))
    jj=which(abs(yc-py)==min(abs(yc-py)))
    if(length(ii)>1 | length(jj)>1){
      ii=ii[1]
      jj=jj[1]
    }
    if(mask[ii,jj,nbands]==0){
      gl.i=ii
      gl.j=jj
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
      ii=which(abs(xc-px)==min(abs(xc-px)))
      jj=which(abs(yc-py)==min(abs(yc-py)))

      if(length(ii)>1 | length(jj)>1){
        ii=ii[1]
        jj=jj[1]
      }

      if(mask[ii,jj,t]==0){
        gl.i=ii
        gl.j=jj
        break
      }else{
        # nothing
      }
    }
    d[t]=sqrt((xc[pd.gl[1]]-xc[gl.i])^2 + (yc[pd.gl[2]]-yc[gl.j])^2)
    if(gl.i<pd.gl[1] | gl.j<pd.gl[2]){d[t]=-d[t]}
  }
  
  nom[i,]=ss.norm[i]*d
  print(i)
}

mean.gl=c()
for(t in 1:length(time)){
  mean.gl[t]=sum(nom[,t],na.rm=T)/sum(ss.norm,na.rm=T)
}
file.out=paste(work.fldr,"/scripts/score-weighting/mean-gl.txt",sep="")
write.table(mean.gl, file.out, quote = F, sep = " ", row.names=F)
  
# score-weighted standard deviation
nom=array(0, dim=c(nr.best,nbands))
for(i in 1:nr.best){
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  file = paste(out.fldr,"/",b.best[i,1]-1,"/yelmo2D.nc", sep="")
  nc = nc_open(file)
  mask=ncvar_get(nc, "f_ice")
  xc=ncvar_get(nc, "xc")
  yc= ncvar_get(nc, "yc")
  time=ncvar_get(nc,"time")
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    nom[i,]=NA
    next
  }
  nc_close(nc)

  # where is the GL at the present?
  for(pt in 1:dim(dat)[1]){
    px=round(dat[pt,1])
    py=round(dat[pt,2])
    # from coord to i,j
    ii=which(abs(xc-px)==min(abs(xc-px)))
    jj=which(abs(yc-py)==min(abs(yc-py)))
    if(length(ii)>1 | length(jj)>1){
      ii=ii[1]
      jj=jj[1]
    }
    if(mask[ii,jj,nbands]==0){
      gl.i=ii
      gl.j=jj
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
      ii=which(abs(xc-px)==min(abs(xc-px)))
      jj=which(abs(yc-py)==min(abs(yc-py)))

      if(length(ii)>1 | length(jj)>1){
        ii=ii[1]
        jj=jj[1]
      }

      if(mask[ii,jj,t]==0){
        gl.i=ii
        gl.j=jj
        break
      }else{
        # nothing
      }
    }
    d[t]=sqrt((xc[pd.gl[1]]-xc[gl.i])^2 + (yc[pd.gl[2]]-yc[gl.j])^2)
    if(gl.i<pd.gl[1] | gl.j<pd.gl[2]){d[t]=-d[t]}
  }

  nom[i,]=ss.norm[i]*(d-mean.gl)^2
  print(i)
}

sd.gl=c()
w=(length(ss.norm)-1)/length(ss.norm)
for(t in 1:length(time)){
  sd.gl[t]=sqrt(sum(nom[,t],na.rm=T)/(w*sum(ss.norm,na.rm=T)))
}
file.out=paste(work.fldr,"/scripts/score-weighting/sd-gl.txt",sep="")
write.table(sd.gl, file.out, quote = F, sep = " ", row.names=F)
