# Ice cores elev

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

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
rmse.g39=rmse.ps100=rmse.g92=c()
for (i in 1:nr){
  
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  
  file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  nc = nc_open(file)
  pd.band = nbands
  if(nc$dim$time$len < pd.band){
    cat("Simulation ",file,"crashed \n")
    rmse.g92[i]=NA
    rmse.ps100[i]=NA
    rmse.g39[i]=NA
    next
  }

  fice=ncvar_get(nc, "f_ice")
  fice.ps100=fice[ps100.i, ps100.j,]
  fice.g92=fice[g92.i,g92.j,]
  fice.g39=fice[g39.i,g39.j,]
  time.y=ncvar_get(nc,"time")

  ti=which(time.y==-14000)
  tf=which(time.y==0)
  
  ##############  PS100 from Syring et al., 2020  ###############
  y.deg.i.ps100=which(fice.ps100<1)
  # take first time with PS100 as ice free
  tt.ps100=which(y.deg.i.ps100>ti & y.deg.i.ps100<tf)[1]
  y.deg.time.ps100=time.y[y.deg.i.ps100[tt.ps100]]
  
  ps100.deg.time=c(-10200,-9600)
  
  if(is.na(y.deg.time.ps100)==T){
    d.ps100 = 10000
  }else if(y.deg.time.ps100 >= -10200 & y.deg.time.ps100 <= -9600){
    d.ps100 = 1
  }else{
    d.ps100 = min(abs(y.deg.time.ps100-ps100.deg.time))
  }
  sq.k.ps100 = d.ps100^2
  n.ps100 = length(d.ps100)

  # RMSE
  rmse.ps100[i]=sqrt(sum(as.vector(sq.k.ps100),na.rm=T)/n.ps100)
  #print(i)
  print(paste("ps100 = ",d.ps100,sep=""))

  ###############   G92 from Davies et al., 2022  #################
  y.deg.i.g92=which(fice.g92<1)
  # take first time with G92 as ice free
  tt.g92=which(y.deg.i.g92>ti & y.deg.i.g92<tf)[1]
  y.deg.time.g92=time.y[y.deg.i.g92[tt.g92]]
  
  g92.deg.time=c(-11200,-10800)
  #g92.deg.time=c(-9600,-7900) # second deglaciation

  if(is.na(y.deg.time.g92)==T){
    d.g92 = 10000
  }else if(y.deg.time.g92 >= -11200 & y.deg.time.g92 <= -10800){
  #}else if(y.deg.time.g92 >= -9600 & y.deg.time.g92 <= -7900){
    d.g92 = 1
  }else{
    d.g92 = min(abs(y.deg.time.g92-g92.deg.time))
  }
  sq.k.g92 = d.g92^2
  n.g92 = length(d.g92)

  # RMSE
  rmse.g92[i]=sqrt(sum(as.vector(sq.k.g92),na.rm=T)/n.g92)
  #print(i)
  print(paste("g92 = ",d.g92,sep=""))


  ###############   G39 from Hansen et al., 2022  #################
  y.deg.i.g39=which(fice.g39<1)
  # take first time with G39 as ice free
  tt.g39=which(y.deg.i.g39>ti & y.deg.i.g39<tf)[1]
  y.deg.time.g39=time.y[y.deg.i.g39[tt.g39]]

  g39.deg.time=c(-13300,-11100)
  #g92.deg.time=c(-9600,-7900) # second deglaciation

  if(is.na(y.deg.time.g39)==T){
    d.g39 = 10000
  }else if(y.deg.time.g39 >= -13300 & y.deg.time.g39 <= -11100){
  #}else if(y.deg.time.g92 >= -9600 & y.deg.time.g92 <= -7900){
    d.g39 = 1
  }else{
    d.g39 = min(abs(y.deg.time.g39-g39.deg.time))
  }
  sq.k.g39 = d.g39^2
  n.g39 = length(d.g39)

  # RMSE
  rmse.g39[i]=sqrt(sum(as.vector(sq.k.g39),na.rm=T)/n.g39)
  #print(i)
  print(paste("g39 = ",d.g39,sep=""))
  nc_close(nc)

  print(i)
}

## write results to a file
file.out=paste(work.fldr,"/rmse-ps100.txt",sep="")
write.table(rmse.ps100, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-g92.txt",sep="")
write.table(rmse.g92, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-g39.txt",sep="")
write.table(rmse.g39, file.out, quote = F, sep = " ", row.names=F)













