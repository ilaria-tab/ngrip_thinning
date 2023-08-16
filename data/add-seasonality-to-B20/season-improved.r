#install.packages("ncdf4",repos='http://cran.us.r-project.org', lib="/mnt/lustre/home/itabone/R/x86_64-redhat-linux-gnu-library/4.0.3")
#libs = "/mnt/lustre/home/itabone/R/x86_64-redhat-linux-gnu-library/4.0.3"
#.libPaths(new = libs)
#.libPaths()
library(ncdf4)
#library(raster)
#library(RColorBrewer)
library(signal)

exp="GRL-8KM_CLIM-RECON-B20-S4T-highP"

inputfile=paste("/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/",exp,".nc",sep="")
outputfile = paste("/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/",exp,"-monthly.nc",sep="")

# seasonality from Buizert et al. 2018
file="/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18.nc"
nc=nc_open(file)
tas.b18=ncvar_get(nc, "tas")
t.b18=ncvar_get(nc, "time")
xc=ncvar_get(nc,"xc")
yc=ncvar_get(nc,"yc")
nc_close(nc)

file="/home/itabone/data/ice_data/Greenland/GRL-8KM/RECON/GRL-8KM_CLIM-RECON-B18-present.nc"
nc=nc_open(file)
tas.pd=ncvar_get(nc, "t2m")
nc_close(nc)

# get TB18(t) from anomalies
nt=length(t.b18)
nm=12
tas.new=tas.b18
tas.new[]=0
for(t in 1:nt){
  for(m in 1:nm){
    temp =tas.b18[,,m,t] + tas.pd[,,m]
    tas.new[,,m,t]=temp
  }
}

# new dT wrt to present (1950)
tas.b18=tas.new
tas.b18[]=0
for(t in 1:nt){
  for(m in 1:nm){
    temp =tas.new[,,m,t] - tas.new[,,m,2202]
    tas.b18[,,m,t]=temp
  }
}
#write.table(tas.b18, file="/home/itabone/work/tasb18.txt", col.names=T, row.names=F)
rm(tas.pd)
rm(temp)
rm(tas.new)

# B20 
#file="/home/itabone/data/ice_data/Greenland/GRL-16KM/RECON/GRL-16KM_CLIM-RECON-B20-meanT-meanP.nc"
nc=nc_open(inputfile)
tas.b20=ncvar_get(nc, "tas")
pr.b20=ncvar_get(nc,"pr")
t.b20=ncvar_get(nc, "time")
x2D=ncvar_get(nc,"x2D")
y2D=ncvar_get(nc,"y2D")
lat2D=ncvar_get(nc,"lat2D")
lon2D=ncvar_get(nc,"lon2D")
xc=ncvar_get(nc,"xc")
yc=ncvar_get(nc,"yc")
area=ncvar_get(nc,"area")
border=ncvar_get(nc,"border")
zs=ncvar_get(nc, "zs")
nc_close(nc)

# Detect monthly variability at each point and filter it
var.mon=array(0, dim=c(length(xc),length(yc),12,2207))
print("Detect monthly variability at each point and filter it")
for(x in 1:length(xc)){
  for(y in 1:length(yc)){
    for(m in 1:12){
      temp=tas.b18[x,y,m,]
      bf <- butter(3, 1/50, type="low")
      var.mon[x,y,m,]=filtfilt(bf, temp)
    }
  }
  print(paste("x =",x))
}

#write.table(var.mon, filename="/home/itabone/work/varmon.txt", col.names=T, row.names=F)
rm(temp)


# Get mean annual temperature from filtered variability signal
tann.b18=array(0, dim=c(length(xc),length(yc),2207))
tann.b18[]=0
for(t18 in 1:length(t.b18)){
  for(m in 1:12){
    tann.b18[,,t18] = tann.b18[,,t18] + var.mon[,,m,t18]
  } 
  tann.b18[,,t18]=tann.b18[,,t18]/12
  print(t18)
}

#Sys.sleep(10)

# Get the residuals between monthly variability and annual mean (for each point, each month and each year)
res=array(0,dim=c(length(xc),length(yc),12,2207))
print("Get the residuals between monthly variability and annual mean")
for(x in 1:length(xc)){
  for(y in 1:length(yc)){
    for(t18 in 1:length(t.b18)){
      for(m in 1:12){
      res[x,y,m,t18]=var.mon[x,y,m,t18]-tann.b18[x,y,t18]
      }
    }
  }
  print(paste("x=",x))
}

#Sys.sleep(10)

# Apply the monthly residuals to B20 (only tas.mon)
tas.mon=array(0, dim=c(length(xc),length(yc),12,401))
pr.mon=array(0, dim=c(length(xc),length(yc),12,401))
print(" Apply the monthly residuals to B20")
for(x in 1:length(xc)){
  for(y in 1:length(yc)){
    for(t in 1:length(t.b20)){
      # which B20 time corresponds to B18 time?
      t.b18.i=which(t.b20[t]==t.b18-5)
      # apply seasonality to B20 data
      for(m in 1:12){
        tas.mon[x,y,m,t]=tas.b20[x,y,t] + res[x,y,m,t.b18.i]
      }
    }
  }
  print(paste("x=",x))
}

#Sys.sleep(10)

# Distribute annual precipitation over months
print("Distribute annual precipitation over months")
for(t in 1:length(t.b20)){
  for(m in 1:12){
    pr.mon[,,m,t] = pr.b20[,,t]

  }
 print(t)
}



################################### new B20 monthly file ####################################################################
#file="/home/itabone/data/ice_data/Greenland/GRL-16KM/RECON/GRL-16KM_CLIM-RECON-B18.nc"
#nc=nc_open(file)
#xc=nc$dim[["xc"]]
#yc=nc$dim[["yc"]]
#month=nc$dim[["month"]]
##time=nc$dim[["time"]]
#nc_close(nc)
print("Write variables to the netcdf file")
# Define dimensions
xc.def=ncdim_def("xc",units="kilometers", vals=xc, unlim=FALSE, create_dimvar=TRUE, longname="xc")
yc.def=ncdim_def("yc",units="kilometers", vals=yc, unlim=FALSE, create_dimvar=TRUE, longname="yc")
mon.def=ncdim_def("month",units="", vals=seq(1,12), unlim=FALSE, create_dimvar=TRUE, longname="month")
time.def=ncdim_def("time",units="years", vals=t.b20, unlim=FALSE, create_dimvar=TRUE, longname="time")

# Define variables
varx2D.def=ncvar_def("x2D",units="", dim=list(xc.def,yc.def), missval=-9999, prec="float")
vary2D.def=ncvar_def("y2D",units="", dim=list(xc.def,yc.def), missval=-9999, prec="float")
varlat2D.def=ncvar_def("lat2D",units="", dim=list(xc.def,yc.def), missval=-9999, prec="float")
varlon2D.def=ncvar_def("lon2D",units="", dim=list(xc.def,yc.def), missval=-9999, prec="float")
vararea.def=ncvar_def("area",units="", dim=list(xc.def,yc.def),missval=-9999, prec="float")
varborder.def=ncvar_def("border",units="", dim=list(xc.def,yc.def),missval=-9999, prec="float")
vartas.def=ncvar_def("tas",units="K", dim=list(xc.def,yc.def,mon.def,time.def), missval=-9999, longname="Surface air temperature at 2m height (anomaly)", prec="float")
varpr.def=ncvar_def("pr",units="",dim=list(xc.def,yc.def,mon.def,time.def), missval=-9999, longname="Precipitation fraction", prec="float")
varzs.def=ncvar_def("zs",units="m above sea level", dim=list(xc.def,yc.def,time.def), missval=-9999, longname="Ice elevation", prec="float")


vars = list(varx2D.def, vary2D.def,varlat2D.def,varlon2D.def,vararea.def,varborder.def,vartas.def,varpr.def,varzs.def)

# Create the file
#utputfile = "/home/itabone/data/ice_data/Greenland/GRL-16KM/RECON/GRL-16KM_CLIM-RECON-B20-meanT-meanP-monthly.nc"
ncout <- nc_create(outputfile, vars,force_v4=TRUE)

# Put variables
ncvar_put(ncout, varx2D.def, x2D)
ncvar_put(ncout, vary2D.def, y2D)
ncvar_put(ncout, varlat2D.def, lat2D)
ncvar_put(ncout, varlon2D.def, lon2D)
ncvar_put(ncout, vararea.def, area)
ncvar_put(ncout, varborder.def, border)
ncvar_put(ncout, vartas.def, tas.mon)
ncvar_put(ncout, varpr.def, pr.mon)
ncvar_put(ncout, varzs.def, zs)

# Put global attributes
ncatt_put(ncout,0,"title","Badgeley et al., 2020, Temperature and Precipitation Reconstructions")

nc_close(ncout)

print("Program finished")
