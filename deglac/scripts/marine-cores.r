# Deglaciation time at inner/outer coasts from sediment cores

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: marine-cores")

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
rmse.g39=rmse.ps100198=rmse.ps100=rmse.g92=c()
for (i in 1:nr){
  
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  
  file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.g39[i]=NA
    rmse.ps100[i]=NA
    rmse.g92[i]=NA
    rmse.ps100198[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.g39[i]=NA
      rmse.ps100[i]=NA
      rmse.g92[i]=NA
      rmse.ps100198[i]=NA
      next
    }
  }
  #nc_close(nc)

  #fice=ncvar_get(nc, "f_ice")
  fice=ncvar_get(nc,"mask_bed")
  fice.ps100=fice[ps100.i, ps100.j,]
  fice.g92=fice[g92.i,g92.j,]
  fice.g39=fice[g39.i,g39.j,]
  fice.ps100198=fice[ps100198.i, ps100198.j,]
  time.y=ncvar_get(nc,"time")

  ti=which(time.y==-15000)
  tf=which(time.y==0)
  
  ##############  PS100 from Syring et al., 2020  ###############
  y.deg.i.ps100=which(fice.ps100>4 | fice.ps100==0)
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
  y.deg.i.g92=which(fice.g92>4 | fice.g92==0)
  # take first time with G92 as ice free
  tt.g92=which(y.deg.i.g92>ti & y.deg.i.g92<tf)[1]
  y.deg.time.g92=time.y[y.deg.i.g92[tt.g92]]
  
  g92.deg.time=c(-13400,-11200)
  #g92.deg.time=c(-11200,-10800) # second deglaciation

  if(is.na(y.deg.time.g92)==T){
    d.g92 = 10000
  }else if(y.deg.time.g92 >= -13400 & y.deg.time.g92 <= -11200){
  #}else if(y.deg.time.g92 >= -11200 & y.deg.time.g92 <= -10800){
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
  y.deg.i.g39=which(fice.g39>4 | fice.g39==0)
  # take first time with G39 as ice free
  tt.g39=which(y.deg.i.g39>ti & y.deg.i.g39<tf)[1]
  y.deg.time.g39=time.y[y.deg.i.g39[tt.g39]]

  g39.deg.time=c(-14000,-13300)

  if(is.na(y.deg.time.g39)==T){
    d.g39 = 10000
  }else if(y.deg.time.g39 >= -14000 & y.deg.time.g39 <= -13300){
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
  #nc_close(nc)


  ###############  PS100/198 from LLoyd et al., 2023  #################
  y.deg.i.ps100198=which(fice.ps100198>4 | fice.ps100198==0)
  # take first time with G39 as ice free
  tt.ps100198=which(y.deg.i.ps100198>ti & y.deg.i.ps100198<tf)[1]
  y.deg.time.ps100198=time.y[y.deg.i.ps100198[tt.ps100198]]

  ps100198.deg.time=c(-10900,-10000)

  if(is.na(y.deg.time.ps100198)==T){
    d.ps100198 = 10000
  }else if(y.deg.time.ps100198 >= -10900 & y.deg.time.ps100198 <= -10000){
    d.ps100198 = 1
  }else{
    d.ps100198 = min(abs(y.deg.time.ps100198-ps100198.deg.time))
  }
  sq.k.ps100198 = d.ps100198^2
  n.ps100198 = length(d.ps100198)

  # RMSE
  rmse.ps100198[i]=sqrt(sum(as.vector(sq.k.ps100198),na.rm=T)/n.ps100198)
  #print(i)
  print(paste("ps100198 = ",d.ps100198,sep=""))
  nc_close(nc)

  print(i)

}


rmse.marine=(rmse.g39+rmse.g92+rmse.ps100+rmse.ps100198)/4
#rmse.marine=rmse.g39

## write results to a file
file.out=paste(work.fldr,"/rmse-marine.txt",sep="")
write.table(rmse.marine, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-ps100.txt",sep="")
write.table(rmse.ps100, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-g92.txt",sep="")
write.table(rmse.g92, file.out, quote = F, sep = " ", row.names=F)
#
file.out=paste(work.fldr,"/rmse-g39.txt",sep="")
write.table(rmse.g39, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-ps100198.txt",sep="")
write.table(rmse.ps100198, file.out, quote = F, sep = " ", row.names=F)












