# Ice cores elev

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)

print("Start: larsen-morains")

out.fldr=paste("/home/titan/gwgk/gwgk005h/models/yelmo/v1.801/yelmox/output/8km/lhs/B18/deglac/B20")
work.fldr=paste("/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20")
nbands=90 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/load-ice-core-data-8km.r")
#source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-marine-cores.r")
source("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ne-cosmogenic-exposures.r")

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
rmse.blaso=rmse.lambland1=rmse.lambland2=rmse.zi=rmse.sondrem=rmse.blochn=rmse.storoen=rmse.bourb=rmse.kapam=c()
for (i in 1:nr){
  
  #file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  
  file = paste(out.fldr,"/",i-1,"/yelmo2D.nc", sep="")

  # has the sim crashed?
  io=file.exists(file)
  if(io==F){
  cat("Simulation ",file,"did not start \n")
    rmse.blaso[i]=NA
    rmse.lambland1[i]=NA
    rmse.lambland2[i]=NA
    rmse.zi[i]=NA
    rmse.sondrem[i]=NA
    rmse.blochn[i]=NA
    rmse.storoen[i]=NA
    rmse.bourb[i]=NA
    rmse.kapam[i]=NA
    next
  }else{
    nc = nc_open(file)
    pd.band = nbands
    if(nc$dim$time$len < pd.band){
      cat("Simulation ",file,"crashed \n")
      rmse.blaso[i]=NA
      rmse.lambland1[i]=NA
      rmse.lambland2[i]=NA
      rmse.zi[i]=NA
      rmse.sondrem[i]=NA
      rmse.blochn[i]=NA
      rmse.storoen[i]=NA
      rmse.bourb[i]=NA
      rmse.kapam[i]=NA
      next
    }
  }
  #nc_close(nc)

  #fice=ncvar_get(nc, "f_ice")
  fice=ncvar_get(nc,"mask_bed")
  fice.blaso=fice[blaso.i, blaso.j,]
  fice.lambland1=fice[lambland1.i, lambland1.j,]
  fice.lambland2=fice[lambland2.i, lambland2.j,]
  fice.zi=fice[zi.i, zi.j,]
  fice.sondrem=fice[sondrem.i, sondrem.j,]
  fice.blochn=fice[blochn.i, blochn.j,]
  fice.storoen=fice[storoen.i,storoen.j,]
  fice.bourb=fice[bourb.i,bourb.j,]
  fice.kapam=fice[kapam.i,kapam.j,]

  time.y=ncvar_get(nc,"time")

  ti=which(time.y==-15000)
  tf=which(time.y==0)
  
  ##############  INNER SITES from Larsen et al. 2018 and Bennike and Weidick 2001 (Blaso, Lambert Land, ZI, Sondre Mellemland, Bloch Nunatakker, Storstrommer, Midgaardsormen   ###############
  
  # Blaso
  y.deg.blaso=which(fice.blaso>5 | fice.blaso<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.blaso=which(y.deg.blaso>ti & y.deg.blaso<tf)[1]
  y.deg.time.blaso=time.y[y.deg.blaso[tt.blaso]]
  
  blaso.deg.time=c(blaso.age.min, blaso.age.max)
  
  if(is.na(y.deg.time.blaso)==T){
    d.blaso = 10000
  }else if(y.deg.time.blaso >=blaso.age.min  & y.deg.time.blaso <= blaso.age.max){
    d.blaso = 1
  }else{
    d.blaso = min(abs(y.deg.time.blaso-blaso.deg.time))
  }
  sq.k.blaso = d.blaso^2
  n.blaso = length(d.blaso)

  # RMSE
  rmse.blaso[i]=sqrt(sum(as.vector(sq.k.blaso),na.rm=T)/n.blaso)
  #print(i)
  print(paste("blaso = ",d.blaso,sep=""))


  ###############   Lambert Land  #################
  y.deg.lambland1=which(fice.lambland1>5 | fice.lambland1<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.lambland1=which(y.deg.lambland1>ti & y.deg.lambland1<tf)[1]
  y.deg.time.lambland1=time.y[y.deg.lambland1[tt.lambland1]]

  lambland1.deg.time=c(lambland1.age.min, lambland1.age.max)

  if(is.na(y.deg.time.lambland1)==T){
    d.lambland1 = 10000
  }else if(y.deg.time.lambland1 >=lambland1.age.min  & y.deg.time.lambland1 <= lambland1.age.max){
    d.lambland1 = 1
  }else{
    d.lambland1 = min(abs(y.deg.time.lambland1-lambland1.deg.time))
  }
  sq.k.lambland1 = d.lambland1^2
  n.lambland1 = length(d.lambland1)

  # RMSE
  rmse.lambland1[i]=sqrt(sum(as.vector(sq.k.lambland1),na.rm=T)/n.lambland1)
  #print(i)
  print(paste("lambland1 = ",d.lambland1,sep=""))


  ############# Lambert Land 2 ####################
  y.deg.lambland2=which(fice.lambland2>5 | fice.lambland2<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.lambland2=which(y.deg.lambland2>ti & y.deg.lambland2<tf)[1]
  y.deg.time.lambland2=time.y[y.deg.lambland2[tt.lambland2]]

  lambland2.deg.time=c(lambland2.age.min, lambland2.age.max)

  if(is.na(y.deg.time.lambland2)==T){
    d.lambland2 = 10000
  }else if(y.deg.time.lambland2 >=lambland2.age.min  & y.deg.time.lambland2 <= lambland2.age.max){
    d.lambland2 = 1
  }else{
    d.lambland2 = min(abs(y.deg.time.lambland2-lambland2.deg.time))
  }
  sq.k.lambland2 = d.lambland2^2
  n.lambland2 = length(d.lambland2)

  # RMSE
  rmse.lambland2[i]=sqrt(sum(as.vector(sq.k.lambland2),na.rm=T)/n.lambland2)
  #print(i)
  print(paste("lambland2 = ",d.lambland2,sep=""))


  ############# Zachariae Istrom ####################
  y.deg.zi=which(fice.zi>5 | fice.zi<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.zi=which(y.deg.zi>ti & y.deg.zi<tf)[1]
  y.deg.time.zi=time.y[y.deg.zi[tt.zi]]

  zi.deg.time=c(zi.age.min, zi.age.max)

  if(is.na(y.deg.time.zi)==T){
    d.zi = 10000
  }else if(y.deg.time.zi >=zi.age.min  & y.deg.time.zi <= zi.age.max){
    d.zi = 1
  }else{
    d.zi = min(abs(y.deg.time.zi-zi.deg.time))
  }
  sq.k.zi= d.zi^2
  n.zi = length(d.zi)

  # RMSE
  rmse.zi[i]=sqrt(sum(as.vector(sq.k.zi),na.rm=T)/n.zi)
  #print(i)
  print(paste("lambland2 = ",d.zi,sep=""))


  ############# Sondre Mellemland  ####################
  y.deg.sondrem=which(fice.sondrem>5 | fice.sondrem<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.sondrem=which(y.deg.sondrem>ti & y.deg.sondrem<tf)[1]
  y.deg.time.sondrem=time.y[y.deg.sondrem[tt.sondrem]]

  sondrem.deg.time=c(sondrem.age.min, sondrem.age.max)

  if(is.na(y.deg.time.sondrem)==T){
    d.sondrem = 10000
  }else if(y.deg.time.sondrem >=sondrem.age.min  & y.deg.time.sondrem <= sondrem.age.max){
    d.sondrem = 1
  }else{
    d.sondrem = min(abs(y.deg.time.sondrem-sondrem.deg.time))
  }
  sq.k.sondrem= d.sondrem^2
  n.sondrem = length(d.sondrem)

  # RMSE
  rmse.sondrem[i]=sqrt(sum(as.vector(sq.k.sondrem),na.rm=T)/n.sondrem)
  #print(i)
  print(paste("sondrem = ",d.sondrem,sep=""))



  ############# Bloch Nunatakker  ####################
  y.deg.blochn=which(fice.blochn>5 | fice.blochn<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.blochn=which(y.deg.blochn>ti & y.deg.blochn<tf)[1]
  y.deg.time.blochn=time.y[y.deg.blochn[tt.blochn]]

  blochn.deg.time=c(blochn.age.min, blochn.age.max)

  if(is.na(y.deg.time.blochn)==T){
    d.blochn = 10000
  }else if(y.deg.time.blochn >=blochn.age.min  & y.deg.time.blochn <= blochn.age.max){
    d.blochn = 1
  }else{
    d.blochn = min(abs(y.deg.time.blochn-blochn.deg.time))
  }
  sq.k.blochn = d.blochn^2
  n.blochn = length(d.blochn)

  # RMSE
  rmse.blochn[i]=sqrt(sum(as.vector(sq.k.blochn),na.rm=T)/n.blochn)
  #print(i)
  print(paste("blochn = ",d.blochn,sep=""))


 ################## OUTER COAST ######################
 ####################################################

  ################# BOURBON OER ######################
  y.deg.bourb=which(fice.bourb>5 | fice.bourb<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.bourb=which(y.deg.bourb>ti & y.deg.bourb<tf)[1]
  y.deg.time.bourb=time.y[y.deg.bourb[tt.bourb]]

  bourb.deg.time=c(bourb.age.min, bourb.age.max)

  if(is.na(y.deg.time.bourb)==T){
    d.bourb = 10000
  }else if(y.deg.time.bourb >=bourb.age.min  & y.deg.time.bourb <= bourb.age.max){
    d.bourb = 1
  }else{
    d.bourb = min(abs(y.deg.time.bourb-bourb.deg.time))
  }
  sq.k.bourb = d.bourb^2
  n.bourb = length(d.bourb)

  # RMSE
  rmse.bourb[i]=sqrt(sum(as.vector(sq.k.bourb),na.rm=T)/n.bourb)
  #print(i)
  print(paste("bourb = ",d.bourb,sep=""))


  ################# Storoen ######################
  y.deg.storoen=which(fice.storoen>5 | fice.storoen<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.storoen=which(y.deg.storoen>ti & y.deg.storoen<tf)[1]
  y.deg.time.storoen=time.y[y.deg.storoen[tt.storoen]]

  storoen.deg.time=c(storoen.age.min, storoen.age.max)

  if(is.na(y.deg.time.storoen)==T){
    d.storoen = 10000
  }else if(y.deg.time.storoen >=storoen.age.min  & y.deg.time.storoen <= storoen.age.max){
    d.storoen = 1
  }else{
    d.storoen = min(abs(y.deg.time.storoen-storoen.deg.time))
  }
  sq.k.storoen = d.storoen^2
  n.storoen = length(d.storoen)

  # RMSE
  rmse.storoen[i]=sqrt(sum(as.vector(sq.k.storoen),na.rm=T)/n.storoen)
  #print(i)
  print(paste("storoen = ",d.storoen,sep=""))


  ################# Kap Amelie ######################
  y.deg.kapam=which(fice.kapam>5 | fice.kapam<=1) #when becomes island, partial, land, or ocean?
  # take first time with PS100 as ice free
  tt.kapam=which(y.deg.kapam>ti & y.deg.kapam<tf)[1]
  y.deg.time.kapam=time.y[y.deg.kapam[tt.kapam]]

  kapam.deg.time=c(kapam.age.min, kapam.age.max)

  if(is.na(y.deg.time.kapam)==T){
    d.kapam = 10000
  }else if(y.deg.time.kapam >=kapam.age.min  & y.deg.time.kapam <= kapam.age.max){
    d.kapam = 1
  }else{
    d.kapam = min(abs(y.deg.time.kapam-kapam.deg.time))
  }
  sq.k.kapam = d.kapam^2
  n.kapam = length(d.kapam)

  # RMSE
  rmse.kapam[i]=sqrt(sum(as.vector(sq.k.kapam),na.rm=T)/n.kapam)
  #print(i)
  print(paste("kapam = ",d.kapam,sep=""))

  nc_close(nc)
  print(i)
}

# MEAN INNER COASTS 
rmse.inner =(rmse.blaso+rmse.lambland1+rmse.lambland2+rmse.zi+rmse.sondrem+rmse.blochn)/6
rmse.outer = (rmse.bourb+rmse.storoen+rmse.kapam)/3

# write results to a file
file.out=paste(work.fldr,"/rmse-inner-coast.txt",sep="")
write.table(rmse.inner, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-outer-coast.txt",sep="")
write.table(rmse.outer, file.out, quote = F, sep = " ", row.names=F)


file.out=paste(work.fldr,"/rmse-kapam.txt",sep="")
write.table(rmse.kapam, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-storoen.txt",sep="")
write.table(rmse.storoen, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse.bourb.txt",sep="")
write.table(rmse.bourb, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-blochn.txt",sep="")
write.table(rmse.blochn, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-sondrem.txt",sep="")
write.table(rmse.sondrem, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-zi.txt",sep="")
write.table(rmse.zi, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-lambland1.txt",sep="")
write.table(rmse.lambland1, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-lambland2.txt",sep="")
write.table(rmse.lambland2, file.out, quote = F, sep = " ", row.names=F)

file.out=paste(work.fldr,"/rmse-blaso.txt",sep="")
write.table(rmse.blaso, file.out, quote = F, sep = " ", row.names=F)
