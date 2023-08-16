# Plots

library(raster)
library(ncdf4)
library(colorRamps)
library(RColorBrewer)
library(plot3D)
library(viridis)

out.fldr=paste("/home/itabone/models/yelmo/v1.23/yelmox/output/8km/test5")
work.fldr=paste("/home/itabone/work/v1.23/lhs/B18/8km/test5")
nbands=42 # how many bands does yelmo2D.nc have?

source(paste(work.fldr,"/scripts/read_LHS.r", sep=""))    # get a, nr, nc, median.par
source(paste(work.fldr,"/scripts/score.r", sep=""))   # get b.best, nr (best)
nr.best=nrow(b.best)

# load surface elevation and climatological data at ice core sites
source("/home/itabone/work/load-ice-core-data-8km.r")

colors=colorRampPalette(rev(c("#3f0d12","#a71d31","#df928e","#edd9a3","#eeffdb","#3ddc97","#52d1dc","#347fc4","#235789")))

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
# Thin rates vs beta

thin.rate=u.negis=taub.negis=taub.egrip=beta.egrip=beta.negis=c()
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  min.i = pd.i

  thin.rate[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000  # m/kyr

  # calc beta averaged over NEGIS
  beta.egrip[i]=log10(mean(beta[egrip.i,egrip.j,max.i:min.i]))
  beta.negis[i]=log10(mean(beta[,,max.i:min.i][basins==9]))
}

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-beta-negis.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(1,7)
ylim=c(-20,0)
plot(xlim, ylim, type="n", xlab="Beta NEGIS log10(Pa yr m-1)", ylab="Mean NGRIP thinning rate (m/kyr)")
points(beta.negis, thin.rate,col="black", bg=mycol,pch=21, cex=1.3)
res=cor.test(beta.negis, thin.rate, method="pearson")
text(2,-17,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(2,-19,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ beta.negis)
beta.new=seq(0,10)
lines(beta.new, beta.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()

##########################################################################################
# Thin rates vs beta (PD)

thin.rate=beta.negis=c()
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  # calc thinning rate
  zs.ngrip=zs[ngrip.i,ngrip.j,]
  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  thin.rate[i]=(zs.ngrip[max.i]-zs.ngrip[pd.i])/(time[max.i]-time[pd.i])*1000  # m/kyr

  # calc beta averaged over NEGIS
  beta.negis[i]=log10(mean(beta[,,pd.i][basins==9]))
}

plot.out=paste(work.fldr,"/Figures/thin-rates-vs-beta-pd.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(3,10)
ylim=c(-20,0)
plot(xlim, ylim, type="n", xlab="Beta NEGIS PD (log10(Pa yr m-1))", ylab="Mean NGRIP thinning rate (m/kyr)")
points(beta.negis, thin.rate,col="black", bg=mycol,pch=21, cex=1.3)
res=cor.test(beta.negis, thin.rate, method="pearson")
text(4,-17,paste("r = ",round(res$estimate[[1]], digits=2),sep=""))
text(4,-19,paste("p = ",signif(res$p.value, digits=4),sep=""),cex=0.7)
# linear regression
lin.reg=lm(thin.rate ~ beta.negis)
beta.new=seq(0,10)
lines(beta.new, beta.new*lin.reg$coefficients[[2]] + lin.reg$coefficients[[1]], col="grey30")
dev.off()


##########################################################################################
# u_NGRIP vs u_NEGIS
plot.out=paste(work.fldr,"/Figures/uNGRIP-vs-uEGRIP.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,80)
ylim=c(0,10)
plot(xlim, ylim, type="n", xlab="U EGRIP (m/yr)", ylab="U NGRIP (m/yr)")
#u.tot=array(0,dim=c(42,2,nr))
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")
  nc_close(nc)

  max.i=which(zs.ngrip==max(zs.ngrip))
  pd.i=length(time)
  u.egrip=u[egrip.i,egrip.j,max.i:pd.i]
  u.ngrip=u[ngrip.i,ngrip.j,max.i:pd.i]
  u.negis=c()
  j=0
  for(t in max.i:pd.i){
   j = j+1 
   u.negis[j]=mean(u[,,t][basins==9])
  }

  points(u.egrip, u.ngrip, pch=21, col="black",bg=mycol[i], cex=1.3)
}
dev.off()

##################################################################################
# uerr_NEGIS vs taub_NEGIS
plot.out=paste(work.fldr,"/Figures/uerrNEGIS-vs-taubNEGIS.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(40000,80000)
ylim=c(0,30)

plot(xlim, ylim, type="n", xlab="Tb (Pa)", ylab="Uerr NEGIS (m/yr)")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, uerr.negis, pch=21,cex=1.3, col="black", bg=mycol[i])
}
dev.off()

###############################################################################
# uerr_NEGIS vs beta_NEGIS
plot.out=paste(work.fldr,"/Figures/uerrNEGIS-vs-betaNEGIS.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,10e6)
ylim=c(0,30)
plot(xlim, ylim, type="n", xlab="Beta (Pa yr/m)", ylab="Uerr NEGIS (m/yr)")
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  beta.negis=mean(beta[,,pd.i][basins==9])
  points(beta.negis, uerr.negis, pch=21,cex=1.3,col="black", bg=mycol[i])
}
dev.off()


################################################################################
# uerr_NEGIS vs q

plot.out=paste(work.fldr,"/Figures/uerrNEGIS-vs-q.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,1)
ylim=c(0,30)
plot(xlim, ylim, type="n", xlab="q", ylab="Uerr NEGIS (m/yr)")
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.negis=mean(taub[,,pd.i][basins==9])
  points(q[i], uerr.negis, pch=21, col="black",bg=mycol[i],cex=1.3)
}
dev.off()

################################################################################
# uerr_NEGIS vs cf_stream

plot.out=paste(work.fldr,"/Figures/uerrNEGIS-vs-cf.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(0,0.5)
ylim=c(0,30)
plot(xlim, ylim, type="n", xlab="cf_NEGIS", ylab="Uerr NEGIS (m/yr)")
for(i in 1:nr.best){
   file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  time=ncvar_get(nc,"time")

  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/length(basins==9))
  taub.negis=mean(taub[,,pd.i][basins==9])
  points(cfngs[i], uerr.negis, pch=21,col="black",bg=mycol[i],cex=1.3)
}
dev.off()




##############################################################################

# Herr_NEGIS vs taub_NEGIS
plot.out=paste(work.fldr,"/Figures/herrNEGIS-vs-taubNEGIS.png", sep="")
png(plot.out, width = 5, height = 5, units = "in", res=100)
#dev.new(width = 5, height = 5, units="in")
xlim=c(20000,100000)
ylim=c(0,50)

plot(xlim, ylim, type="n", xlab="Tb NEGIS (Pa)", ylab="Herr NEGIS (m)")
for(i in 1:nr.best){
  file = paste(out.fldr,"/itmb.",itmb[i],".itmc.",itmc[i],".btq.",q[i],".cbz0.",z0[i],".cfngs.",cfngs[i],".nffdlt.",neff[i],"/yelmo2D.nc", sep="")
  nc=nc_open(file)
  basins=ncvar_get(nc,"basins")
  beta=ncvar_get(nc,"beta")
  taub=ncvar_get(nc,"taub")
  uerr=ncvar_get(nc,"uxy_s_pd_err")
  u=ncvar_get(nc,"uxy_s")
  zs=ncvar_get(nc,"z_srf")
  herr=ncvar_get(nc,"H_ice_pd_err")
  time=ncvar_get(nc,"time")

  uerr.negis=sqrt(sum((uerr[,,pd.i][basins==9])^2)/(length(basins==9)))
  herr.negis=sqrt(sum((herr[,,pd.i][basins==9])^2)/(length(basins==9)))
  taub.negis=mean(taub[,,max.i:min.i][basins==9])
  points(taub.negis, herr.negis, pch=21,col="black",bg=mycol[i],cex=1.3)
}
dev.off()


####################################################################################################################

# zsNGRIP, thin rate (derivative), gl retreat vs time

#for(run in 1:1){
#   file = paste(out.fldr,"/itmb.",itmb[run],".itmc.",itmc[run],".btq.",q[run],".cbz0.",z0[run],".cfngs.",cfngs[run],".nffdlt.",neff[run],"/yelmo2D.nc", sep="")
#  nc=nc_open(file)
#  mask=ncvar_get(nc, "mask_bed")
#  zs=ncvar_get(nc,"z_srf")
#  u=ncvar_get(nc,"uxy_s")
#  xc=ncvar_get(nc, "xc")
#  yc= ncvar_get(nc, "yc")
#  time=ncvar_get(nc,"time")
#  dtann=ncvar_get(nc,"dTa_ann")
#  dtann.neoff=dtann[85,150,]
#    
#  # u at EGRIP
#  u.egrip=u[egrip.i,egrip.j,]
#  acc.egrip=c()
#  for(t in 11:length(time)){
#    acc.egrip[t]=(u.egrip[t]-u.egrip[t-1])/(time[t]-time[t-1])*1000
#  }
#
#  # where is the GL at the present?
#  for(pt in 1:dim(dat)[1]){
#    px=round(dat[pt,1])
#    py=round(dat[pt,2])
#    # from coord to i,j
#    i=which(abs(xc-px)==min(abs(xc-px)))
#    j=which(abs(yc-py)==min(abs(yc-py)))
#    if(length(i)>1 | length(j)>1){
#      i=i[1]
#      j=j[1]
#    }
#    if(mask[i,j,42]==4){
#      gl.i=i
#      gl.j=j
#      break
#    }else{
#    }
#  }
#  pd.gl=c(gl.i,gl.j)
#
#  # calc distance from pd.gl
#  d=c()
#  for(t in 1:length(time)){
#    for(pt in 1:dim(dat)[1]){
#      px=round(dat[pt,1])
#      py=round(dat[pt,2])
#      # from coord to i,j
#      i=which(abs(xc-px)==min(abs(xc-px)))
#      j=which(abs(yc-py)==min(abs(yc-py)))
#      if(length(i)>1 | length(j)>1){
#        i=i[1]
#        j=j[1]
#      }
#      if(mask[i,j,t]==4){
#        gl.i=i
#        gl.j=j
#        break
#      }else{
#      }
#    }
#    d[t]=sqrt((xc[pd.gl[1]]-xc[gl.i])^2 + (yc[pd.gl[2]]-yc[gl.j])^2)
#    if(gl.i<pd.gl[1] | gl.j<pd.gl[2]){d[t]=-d[t]}
#  }
#  # retreat rates
#  pd.i=length(time)
#  retr=c()
#  for(t in 11:length(time)){
#    retr[t]=(d[t]-d[t-1])/(time[t]-time[t-1]) *1000
#  }
#
#  # retreat acceleration
#  retr.acc=c()
#  for(t in 11:length(time)){
#    retr.acc[t]=(retr[t]-retr[t-1])/(time[t]-time[t-1]) *1000
#  }
#
#  # thinning rates
#  zs.ngrip=zs[ngrip.i,ngrip.j,]
#  dzs.ngrip=zs.ngrip-zs.ngrip[pd.i]
#  max.i=which(zs.ngrip==max(zs.ngrip))
#  rate=c()
#  for(t in 11:length(time)){
#    rate[t]=(zs.ngrip[t]-zs.ngrip[t-1])/(time[t]-time[t-1])*1000
#  }
#  
#  #plot.out=paste(work.fldr,"/Figures/1d-",run,".png", sep="")
#  #png(plot.out, width = 6, height = 6, units = "in", res=100)
#  dev.new(width=6,height=6,units="in")
#  par(mfrow=c(2,1), mar=c(0,3,0,3), oma=c(4,1,1,1))
#  
#  # NGRIP thin rate
#  at.x=seq(-16000,0,by=2000)
#  at.y=seq(-50,50, by=25)
#  xlim=c(-16000,0)
#  ylim=c(-50,50)
#  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
#  lines(time, -rate,col="#342E37", lwd=3)
#  axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
#  mtext("NGRIP thinnig rate (m/kyr)",side=2,line=2.7,cex=1)
#  
#  # EGRIP U
#  at.y=seq(-10,50, by=20)
#  labels.y=seq(0,60,by=20)
#  lines(time, u.egrip-10, col="#6CC25B", lwd=3)
#  axis(4,at=at.y,labels=labels.y,col="#6CC25B",col.axis="#6CC25B",las=1, cex.axis=1)
#  mtext("          U  EGRIP (m/yr)",side=4,line=2.7,cex=1)
#  
#  # GL distance
#  at.y=seq(-100,350, by=50)
#  ylim=c(-100,350)
#  plot(xlim, ylim, type="n", ann=F, axes=F, xaxt="n", cex.axis=1)
#  lines(time, d, col="#E93835", lwd=3)
#  axis(2,at=at.y,labels=at.y,col="black",col.axis="black",las=1, cex.axis=1)
#  mtext("NEGIS GL retreat (km)",side=2,line=2.7,cex=1)
#  
#  # dTann at NE 
#  at.y=seq(-40,240, by=40)
#  labels.y=seq(-30,5,by=5)
#  lines(time, dtann.neoff*8+200, col="#1098F7",lwd=3)
#  axis(4,at=at.y,labels=labels.y,col="black",col.axis="black",las=1, cex.axis=1)
#  mtext("NE dTann (°C)",side=4,line=2.7,cex=1)
#  axis(1,at=at.x, labels=at.x, col="black", col.axis="black", las=1, cex.axis=1)
#  mtext("Time (yr)",side=1,line=2.5,cex=1)
#  #dev.off()
#}
#
#
