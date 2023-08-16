# ice core locations (yelmo2D index)

fun.file=("/home/titan/gwgk/gwgk005h/work/ngrip/scripts/fun_findij-4km.r")
source(fun.file)

# OUTER COAST SITES 
# Bourbon Oer (Larsen et al., 2018)
lat=78.620
lon=-18.402
bourb.i=findij(lat,lon)[1]
bourb.j=findij(lat,lon)[2]
bourb.x=findij(lat,lon)[3]
bourb.y=findij(lat,lon)[4]

# age (ka)
age=c(-10.79,-11.44,-11.03)*1000
err=c(0.39,0.35,0.37)*1000
bourb.age.min=min(age-err)
bourb.age.max=max(age+err)

# Storoen (Larsen et al., 2018)
lat=78.066
lon=-19.091
storoen.i=findij(lat,lon)[1]
storoen.j=findij(lat,lon)[2]
storoen.x=findij(lat,lon)[3]
storoen.y=findij(lat,lon)[4]

age=c(-11.3,-11.6,-11.26)*1000
err=c(0.33,0.85,0.42)*1000
storoen.age.min=min(age-err)
storoen.age.max=max(age+err)

# Kap Amelie (Larsen et al., 2018)
lat=77.544
lon=-19.131
kapam.i=findij(lat,lon)[1]
kapam.j=findij(lat,lon)[2]
kapam.x=findij(lat,lon)[3]
kapam.y=findij(lat,lon)[4]

age=c(-12.97,-10.34,-14.2)*1000
err=c(0.73,1.04,0.46)*1000
kapam.age.min=min(age-err)
kapam.age.max=max(age+err)


# INNER (PRESENT-DAY) SITES
# Blaso (Larsen et al., 2018)
lat=79.636
lon=-23.104
blaso.i=findij(lat,lon)[1]
blaso.j=findij(lat,lon)[2]
blaso.x=findij(lat,lon)[3]
blaso.y=findij(lat,lon)[4]

age=c(-10.22,-7.86,-12.11)*1000
err=c(0.42,2.69,1.47)*1000
blaso.age.min=min(age-err)
blaso.age.max=max(age+err)

# Lambert Land 1 (Larsen et al., 2018)
lat=79.143
lon=-21.41
lambland1.i=findij(lat,lon)[1]
lambland1.j=findij(lat,lon)[2]
lambland1.x=findij(lat,lon)[3]
lambland1.y=findij(lat,lon)[4]

age=c(-9.27,-9.08,-8.96)*1000
err=c(0.46,0.56,0.65)*1000

lambland1.age.min=min(age-err)
lambland1.age.max=max(age+err)


# Lambert Land 2 (Larsen et al., 2018)
lat= 79.1
lon=-20.9
lambland2.i=findij(lat,lon)[1]
lambland2.j=findij(lat,lon)[2]
lambland2.x=findij(lat,lon)[3]
lambland2.y=findij(lat,lon)[4]

age=c(-8.99,-8.65)*1000
err=c(0.49,0.54)*1000

lambland2.age.min=min(age-err)
lambland2.age.max=max(age+err)


# Zachariae Istrom (Larsen et al., 2018)
lat=78.782
lon=-20.777
zi.i=findij(lat,lon)[1]
zi.j=findij(lat,lon)[2]
zi.x=findij(lat,lon)[3]
zi.y=findij(lat,lon)[4]

age=c(-8.83,-8.55,-9.86)*1000
err=c(0.5,0.41,0.48)*1000
zi.age.min=min(age-err)
zi.age.max=max(age+err)


# Sondre Mellemland (Larsen et al., 2018)
lat=78.067
lon=-21.505
sondrem.i=findij(lat,lon)[1]
sondrem.j=findij(lat,lon)[2]
sondrem.x=findij(lat,lon)[3]
sondrem.y=findij(lat,lon)[4]

age=c(-10.82,-9.03,-9.36)*1000
err=c(0.67,0.45,0.5)*1000

sondrem.age.min=min(age-err)
sondrem.age.max=max(age+err)

# Bloch Nunatakker (Larsen et al., 2018)
lat=79.6
lon=-19.524
blochn.i=findij(lat,lon)[1]
blochn.j=findij(lat,lon)[2]
blochn.x=findij(lat,lon)[3]
blochn.y=findij(lat,lon)[4]

age=c(-9.07,-8.59,-9.21)*1000
err=c(0.35,0.32,0.62)*1000

blochn.age.min=min(age-err)
blochn.age.max=max(age+err)


#############  Bennike and Weidick 2001
## Blaso
#
#degrees.to.decimals<-function(degrees=45,minutes=30)
#{
#  if(!is.numeric(minutes)) stop("Please enter a numeric value for minutes!\n")
#  if(!is.numeric(degrees)) stop("Please enter a numeric value for degrees!\n")
#  decimal=minutes/60
#  result=degrees+decimal
#  result
#}
#
#lat=c(degrees.to.decimals(79,34.6), degrees.to.decimals(79,34.9), degrees.to.decimals(79,37.6),
#      degrees.to.decimals(79,36.8), degrees.to.decimals(79,37.6), degrees.to.decimals(79,36.8),
#      degrees.to.decimals(79,37.2), degrees.to.decimals(79,40), degrees.to.decimals(79,37.6),
#      degrees.to.decimals(79,34.6), degrees.to.decimals(79,34.3), degrees.to.decimals(79,32.9),
#      degrees.to.decimals(79,36.8))
#lon=c(degrees.to.decimals(22,37.3), degrees.to.decimals(22,18.2), degrees.to.decimals(21,55.5), 
#      degrees.to.decimals(22,30.9), degrees.to.decimals(22,23.0), degrees.to.decimals(22,30.9), 
#      degrees.to.decimals(22,35.0), degrees.to.decimals(22,30), degrees.to.decimals(22,23),
#      degrees.to.decimals(22,37.3), degrees.to.decimals(22,23.3), degrees.to.decimals(22,37.4),
#      degrees.to.decimals(22,30.9))
#lat=mean(lat)
#lon=mean(lon)
#
#in.bla.bw.i=findij(lat,lon)[1]
#in.bla.bw.j=findij(lat,lon)[2]
#in.bla.bw.x=findij(lat,lon)[3]
#in.bla.bw.y=findij(lat,lon)[4]
#
#age=c(-6940,-6755,-6080,-6040,-6035,-5995,-5955,-5440,-5080,-5585,-5190,-4645,-4590)
#err=c(75,55,55,55,55,55,55,55,55,65,70,55,45)
#in.bla.bw.min=min(age-err)
#in.bla.bw.max=max(age+err)
#
## Midgaardsormen 
#lat=c(degrees.to.decimals(79,35.2), degrees.to.decimals(79,39.7), degrees.to.decimals(79,32.3),
#      degrees.to.decimals(79,39.7), degrees.to.decimals(79,39.8), degrees.to.decimals(79,39.7),
#      degrees.to.decimals(79,38.9), degrees.to.decimals(79,38.9), degrees.to.decimals(79,38.9),
#      degrees.to.decimals(79,38.9), degrees.to.decimals(79,38.9), degrees.to.decimals(79,36.8),
#      degrees.to.decimals(79,39.5), degrees.to.decimals(79,39.5), degrees.to.decimals(79,34.4), 
#      degrees.to.decimals(79,36.6), degrees.to.decimals(79,39.5), degrees.to.decimals(79,35.2),
#      degrees.to.decimals(79,39.5), degrees.to.decimals(79,39.0), degrees.to.decimals(79,39.6), 
#      degrees.to.decimals(79,38.9))
#lon=c(degrees.to.decimals(21,40.9), degrees.to.decimals(21,5.2), degrees.to.decimals(22,26.6),
#      degrees.to.decimals(21,5.2), degrees.to.decimals(21,2.6), degrees.to.decimals(21,2.8),
#      degrees.to.decimals(21,9.0), degrees.to.decimals(21,9.2), degrees.to.decimals(21,9.9),
#      degrees.to.decimals(21,9.0), degrees.to.decimals(21,9.0), degrees.to.decimals(22,30.9),
#      degrees.to.decimals(21,8.8), degrees.to.decimals(21,1.8), degrees.to.decimals(22,23.3), 
#      degrees.to.decimals(22,18.3), degrees.to.decimals(21,1.8), degrees.to.decimals(21,40.9), 
#      degrees.to.decimals(21,1.8), degrees.to.decimals(21,7.2), degrees.to.decimals(21,3.0), 
#      degrees.to.decimals(21,8.1))
#lat=mean(lat)
#lon=mean(lon)
#
#in.mid.bw.i=findij(lat,lon)[1]
#in.mid.bw.j=findij(lat,lon)[2]
#in.mid.bw.x=findij(lat,lon)[3]
#in.mid.bw.y=findij(lat,lon)[4]
#
#age=c(7480,7045,6805,6160,6665,6535,6245,6310,6125,6115,6065,6050,5985,5945,5830,5775,5695,5660,5405,4730,5135,4965)
#err=c(170,75,65,75,70,80,80,75,80,75,80,105,65,70,50,60,80,90,85,95,70,55)
#in.mid.bw.min=min(age-err)
#in.mid.bw.max=max(age+err)
#
## Sondre Mellemland
#lat=degrees.to.decimals(78,4.6)
#lon=c(degrees.to.decimals(21,38), degrees.to.decimals(21,42), degrees.to.decimals(21,38),
#      degrees.to.decimals(21,58), degrees.to.decimals(21,58), degrees.to.decimals(21,42))
#lat=mean(lat)
#lon=mean(lon)
#
#in.sm.bw.i=findij(lat,lon)[1]
#in.sm.bw.j=findij(lat,lon)[2]
#in.sm.bw.x=findij(lat,lon)[3]
#in.sm.bw.y=findij(lat,lon)[4]
#
#age=c(5875,5750,5630,5550,5500,5350)
#err=c(60,60,60,55,55,55)
#in.sm.bw.min=min(age-err)
#in.sm.bw.max=max(age+err)
#
## Storstrommen
#lat=c(degrees.to.decimals(77,10), degrees.to.decimals(77,5.1),degrees.to.decimals(77,9.8), degrees.to.decimals(77,10),degrees.to.decimals(77,10),degrees.to.decimals(77,11),degrees.to.decimals(77,11),degrees.to.decimals(77,10),degrees.to.decimals(77,9.9),degrees.to.decimals(77,10),degrees.to.decimals(77,10),degrees.to.decimals(77,10),degrees.to.decimals(77,10))
#lon=c(degrees.to.decimals(21,58),degrees.to.decimals(21,55.2),degrees.to.decimals(21,58.7), degrees.to.decimals(21,58), degrees.to.decimals(21,55),degrees.to.decimals(21,57), degrees.to.decimals(21,57), degrees.to.decimals(22,0),degrees.to.decimals(21,58.7),degrees.to.decimals(22,0),degrees.to.decimals(22,0),degrees.to.decimals(22,0),degrees.to.decimals(22,0))
#lat=mean(lat)
#lon=mean(lon)
#
#in.stor.bw.i=findij(lat,lon)[1]
#in.stor.bw.j=findij(lat,lon)[2]
#in.stor.bw.x=findij(lat,lon)[3]
#in.stor.bw.y=findij(lat,lon)[4]
#
#age=c(1815,3230,3630,3725,4180,4840,4910,5030,5180)
#err=c(55,85,90,60,60,90,85,75,95)
#in.stor.bw.min=min(age-err)
#in.stor.bw.max=max(age+err)
#
