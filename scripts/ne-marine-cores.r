# ice core locations (yelmo2D index)

fun.file=("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-8km.r")
source(fun.file)

# core PS100/270-1 (https://doi.org/10.1594/PANGAEA.921185)
lat=79.497170
lon=-18.139670
ps100.i=findij(lat,lon)[1]
ps100.j=findij(lat,lon)[2]
ps100.x=findij(lat,lon)[3]
ps100.y=findij(lat,lon)[4]

# core DA17-NG-ST-08-092G (https://doi.pangaea.de/10.1594/PANGAEA.945987)
lat=78.500900
lon=-17.278517
g92.i=findij(lat,lon)[1]
g92.j=findij(lat,lon)[2]
g92.x=findij(lat,lon)[3]
g92.y=findij(lat,lon)[4]

