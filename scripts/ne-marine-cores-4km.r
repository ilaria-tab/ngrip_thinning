# ice core locations (yelmo2D index)

fun.file=("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-4km.r")
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

# core DA17-NG-ST03-039G (https://doi.pangaea.de/10.1594/PANGAEA.943880) (Hansen et al., 2022, https://doi.org/10.1016/j.quascirev.2022.107704)
lat=80.037167
lon=-8.923183
g39.i=findij(lat,lon)[1]
g39.j=findij(lat,lon)[2]
g39.x=findij(lat,lon)[3]
g39.y=findij(lat,lon)[4]

# Core in Fram Strait (Werner et al., 2016)
lat=79.161001
lon=5.337833
fs.i=findij(lat,lon)[1]
fs.j=findij(lat,lon)[2]
fs.x=findij(lat,lon)[3]
fs.y=findij(lat,lon)[4]

# Core PS10/198 (Lloyd et al., 2023)
lat=79.191233
lon=-17.107167
ps100198.i=findij(lat,lon)[1]
ps100198.j=findij(lat,lon)[2]
ps100198.x=findij(lat,lon)[3]
ps100198.y=findij(lat,lon)[4]
