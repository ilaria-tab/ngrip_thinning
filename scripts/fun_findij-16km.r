library(ncdf4)
library(raster)
library(fields)

# GRL-16KM
findij=function(lat, lon)
{
  file="/home/titan/gwgk/gwgk005h/data/ice_data/Greenland/GRL-16KM/GRL-16KM_TOPO-M17.nc"
  nc=nc_open(file)
  xc=ncvar_get(nc, "xc")
  yc=ncvar_get(nc, "yc")
  lat2D=ncvar_get(nc, "lat2D")
  lon2D=ncvar_get(nc, "lon2D")
  
  delta=1e5
  
  for(i in 1:length(xc)){
    for(j in 1:length(yc)){
      if( (sqrt( (lat2D[i,j]-lat)^2 + (lon2D[i,j]-lon)^2 ) < delta) ){
        delta = sqrt( (lat2D[i,j]-lat)^2 + (lon2D[i,j]-lon)^2 )
        core.i=i
        core.j=j
      }
    }
  }
  
  #coordinates
  core.x=xc[core.i]
  core.y=yc[core.j]
  
  ij=c(core.i,core.j,core.x,core.y)
  return(ij)
}

  
