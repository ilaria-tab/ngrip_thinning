# ice core locations (yelmo2D index)

fun.file=("/mnt/lustre/scratch/itabone/ilaria/snowball/scripts/fun_findij.r")
source(fun.file)

# NGRIP
lat=75.10 
lon=-42.32
ngrip.i=findij(lat,lon)[1]
ngrip.j=findij(lat,lon)[2]
ngrip.x=findij(lat,lon)[3]
ngrip.y=findij(lat,lon)[4]
  
# GRIP
lat=72.58
lon=-37.64
grip.i=findij(lat,lon)[1]
grip.j=findij(lat,lon)[2]
grip.x=findij(lat,lon)[3]
grip.y=findij(lat,lon)[4]

# Dye3
lat=65.18
lon=-43.82
dye3.i=findij(lat,lon)[1]
dye3.j=findij(lat,lon)[2]
dye3.x=findij(lat,lon)[3]
dye3.y=findij(lat,lon)[4]

# Camp Century
lat=77.17
lon=-61.13
campcen.i=findij(lat,lon)[1]
campcen.j=findij(lat,lon)[2]
campcen.x=findij(lat,lon)[3]
campcen.y=findij(lat,lon)[4]

# NEEM
lat=77.44
lon=-51.07
neem.i=findij(lat,lon)[1]
neem.j=findij(lat,lon)[2]
neem.x=findij(lat,lon)[3]
neem.y=findij(lat,lon)[4]


# GISP2
lat=72.970000
lon=-38.8000
gisp2.i=findij(lat,lon)[1]
gisp2.j=findij(lat,lon)[2]
gisp2.x=findij(lat,lon)[3]
gisp2.y=findij(lat,lon)[4]

# Agassiz
lat=80.7  
lon=-73.1
aga.i=findij(lat,lon)[1]
aga.j=findij(lat,lon)[2]
aga.x=findij(lat,lon)[3]
aga.y=findij(lat,lon)[4]
  
# Renland
lat=71.27
lon=-26.73
ren.i=findij(lat,lon)[1]
ren.j=findij(lat,lon)[2]
ren.x=findij(lat,lon)[3]
ren.y=findij(lat,lon)[4]

# EGRIP
lat=75.630000
lon=-35.99
egrip.i=findij(lat,lon)[1]
egrip.j=findij(lat,lon)[2]
egrip.x=findij(lat,lon)[3]
egrip.y=findij(lat,lon)[4]
