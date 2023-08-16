# NG glacier at 70km from present ice front must never be ice free

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: ng-icefree")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-marine-cores.r")

alpha.col <- function(col,percent=50)
{
  c <- rep("#000000",length(col))
  for (i in 1:length(col)) {
    cc <- col2rgb(col[i])/255
    c[i] <- rgb(t(cc),alpha=percent/100)
  }

  return(c)
}

# Step 2: Calulate misfits
rmse.ng=c()
for (i in 1:nr){
  
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  
  file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.ng[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.ng[i]=NA
      next
    }
  }
  #nc_close(nc)

  #fice=ncvar_get(nc, "f_ice")
  fice=ncvar_get(nc,"mask_bed")
  fice.ng=fice[143, 290,]
  time.y=ncvar_get(nc,"time")

  ti=which(time.y==-15000)
  tf=which(time.y==0)
  
  ##############  NG glacier at 70km from present ice front must never be ice free  ###############
  y.ng=which(fice.ng>5 | fice.ng<=1) # if ocean, land,island or partial
  # take first time with PS100 as ice free
  #tt.ng=which(y.ng>ti & y.ng<tf)[1]
  #y.deg.time.ng=time.y[y.ng[tt.ng]]
  
  deg.time=time.y[y.ng]

  if(any(deg.time >= -15000 & deg.time <= 0)){
    # it deglaciated but it shouldn't 
    d.ng = 1000
  }else{
    # it didn't deglaciated
    d.ng = 1
  }
  
  sq.k.ng = d.ng^2
  n.ng = length(d.ng)

  # RMSE
  rmse.ng[i]=sqrt(sum(as.vector(sq.k.ng),na.rm=T)/n.ng)
  #print(i)
  print(paste("ng = ",d.ng,sep=""))

  nc_close(nc)

  print(i)
}

## write results to a file
file.out=paste(work.fldr,"/rmse-ng.txt",sep="")
write.table(rmse.ng, file.out, quote = F, sep = " ", row.names=F)

#file.out=paste(work.fldr,"/rmse-ps100.txt",sep="")
#write.table(rmse.ps100, file.out, quote = F, sep = " ", row.names=F)

#file.out=paste(work.fldr,"/rmse-g92.txt",sep="")
#write.table(rmse.g92, file.out, quote = F, sep = " ", row.names=F)
#
#file.out=paste(work.fldr,"/rmse-g39.txt",sep="")
#write.table(rmse.g39, file.out, quote = F, sep = " ", row.names=F)













