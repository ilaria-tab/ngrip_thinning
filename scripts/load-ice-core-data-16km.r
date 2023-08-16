library(fields)
library(ncdf4)
library(raster)
library(RColorBrewer)

#########################################################################################################################################################################
# This file loads:
# core locations: ngrip.i, ngrip.j, grip.i, grip.j, dye3.i, dye3.j, campcen.i, campcen.j, neem.i, neem.j, gisp2.i, gisp2.j agas.i, agas.j, ren.i, ren.j, egrip.i, egrip.j
# elevation changes data: Vinther 2009: ds.grip.v09, ds.ngrip.v09, ds.campcen.v09, ds.dye3.v09 (t.ds.v09)
#                         Lecavalier 2013: ds.grip1.l13, ds.grip2.l13, ds.grip.max.l13, ds.grip.min.l13, 
#                                          ds.ngrip1.l13, ds.ngrip2.l13, ds.ngrip.max.l13, ds.ngrip.min.l13,
#                                          ds.campcen1.l13, ds.campcen2.l13, ds.campcen.max.l13, ds.campcen.min.l13, 
#                                          ds.dye31.l13, ds.dye32.l13, ds.dye3.min.l13, ds.dye3.max.l13 (t.ds.l13) 
# climatologic data:      Kindler 2014: acc.ngrip.k14, tann.ngrip.k14 (t.k14)
#                         Cuffey & Clow 1999: acc.gisp2.cc99 (t.cc99)
#                         Alley 2000: acc.gisp2.a00 (t.acc.a00), tann.gisp2.a00 (t.tann.a00)
#                         Rassmussen 2013: acc.neem.mean.r13, acc.neem.min.r13, acc.neem.max.r13 (t.r13)
#                         Badgeley 2020: acc.grip.b20, (t.grip.b20), acc.ngrip.b20 (t.ngrip.b20), acc.dye3.b20 (t.dye3.b20)
#                         Gkinis 2014: acc.ngrip.g14 (t.g14)
#                         Buizert 2018: tann.neem.b18, tann.ngrip.b18, tann.gisp2.b18, tann.dye3.b18, tann.ren.b18, tann.aga.b18, tann.campcen.b18, tann.egrip.b18 (t.b18)
#                                       tsum.neem.b18, tsum.ngrip.b18, tsum.gisp2.b18, tsum.dye3.b18, tsum.ren.b18, tsum.aga.b18, tsum.campcen.b18, tsum.egrip.b18 (t.b18)
###########################################################################################################################################################################

fldr="/home/titan/gwgk/gwgk005h/work/ngrip/data"

# ice core sites
file.cores="/home/titan/gwgk/gwgk005h/work/ngrip/scripts/ice-core-16km.txt"
txt=read.table(file.cores,header=F, sep="", skip=3)
i=txt$V2
j=txt$V3
ngrip.i=i[1]
ngrip.j=j[1]
grip.i=i[2]
grip.j=j[2]
dye3.i=i[3]
dye3.j=j[3]
campcen.i=i[4]
campcen.j=j[4]
neem.i=i[5]
neem.j=j[5]
gisp2.i=i[6]
gisp2.j=j[6]
agas.i=i[7]
agas.j=j[7]
ren.i=i[8]
ren.j=j[8]
egrip.i=i[9]
egrip.j=j[9]

# elevation data from Vinther 2009
filename = paste(fldr,"/vinther2009/Elevation-change_Vinther_2009.txt",sep="")
txt=read.table(filename,header=F, sep="", skip=15)
vinth=as.matrix(txt)
colnames(vinth)=c("yrBP","CC","NGRIP","GRIP","DYE3")

# take elevation change every 100 yr (from PD to 11.7 kyr BP)
j=4
s.change=c()
t=c()
t[1]=0
s.change[1]=0
for(i in 1:117){
  s=vinth[j,3]
  time=vinth[j,1]
  s.change[i+1]=s
  t[i+1]=time
  j=j+5
}
ds.ngrip.v09=rev(s.change)
s.ngrip.v09=ds.ngrip.v09 + 2920

j=4
s.change=c()
t=c()
t[1]=0
s.change[1]=0
for(i in 1:117){
  s=vinth[j,4]
  time=vinth[j,1]
  s.change[i+1]=s
  t[i+1]=time
  j=j+5
}
ds.grip.v09=rev(s.change)
s.grip.v09=ds.grip.v09 + 3230

j=4
s.change=c()
t=c()
t[1]=0
s.change[1]=0
for(i in 1:117){
  s=vinth[j,5]
  time=vinth[j,1]
  s.change[i+1]=s
  t[i+1]=time
  j=j+5
}
ds.dye3.v09=rev(s.change)
s.dye3.v09=ds.dye3.v09 + 2490

j=4
s.change=c()
t=c()
t[1]=0
s.change[1]=0
for(i in 1:117){
  s=vinth[j,2]
  time=vinth[j,1]
  s.change[i+1]=s
  t[i+1]=time
  j=j+5
}
ds.campcen.v09=rev(s.change)  # m
s.campcen.v09=ds.campcen.v09 + 1890

t.ds.v09=-rev(t)  # yr


###### Lecavalier et al., 2013:  NRGIP, GRIP, Dye3, Campcentury Elevation change (8-0 kyr) ########
file=paste(fldr,"/lecavalier2013/Lecavalier_etal_2013_elevation.csv",sep="")
data=read.csv(file, skip=6, sep=";", header=F)
time=rev(data$V1)  # kyr BP
time= time*1000 # yr BP

t.ds.l13=time

# elevation changes for site A and B + min and max values
ds.grip1.l13=rev(data$V2) # m
ds.grip2.l13=rev(data$V3)
ds.grip.max.l13=rev(data$V4)
ds.grip.min.l13=rev(data$V5)

ds.ngrip1.l13=rev(data$V7) # m
ds.ngrip2.l13=rev(data$V8)
ds.ngrip.max.l13=rev(data$V9)
ds.ngrip.min.l13=rev(data$V10)

ds.dye31.l13=rev(data$V12) # m
ds.dye32.l13=rev(data$V13)
ds.dye3.max.l13=rev(data$V14)
ds.dye3.min.l13=rev(data$V15)

ds.campcen1.l13=rev(data$V17) # m
ds.campcen2.l13=rev(data$V18)
ds.campcen.max.l13=rev(data$V19)
ds.campcen.min.l13=rev(data$V20)

# absolute elevations
s.grip1.l13=ds.grip1.l13 + 3230
s.grip2.l13=ds.grip2.l13 + 3230
s.grip.max.l13=ds.grip.max.l13 + 3230
s.grip.min.l13=ds.grip.min.l13 + 3230

s.ngrip1.l13=ds.ngrip1.l13 + 2920
s.ngrip2.l13=ds.ngrip2.l13 + 2920
s.ngrip.max.l13=ds.ngrip.max.l13 + 2920
s.ngrip.min.l13=ds.ngrip.min.l13 + 2920

s.dye31.l13=ds.dye31.l13 + 2490
s.dye32.l13=ds.dye32.l13 + 2490
s.dye3.max.l13=ds.dye3.max.l13 + 2490
s.dye3.min.l13=ds.dye3.min.l13 + 2490

s.campcen1.l13=ds.campcen1.l13 + 1890
s.campcen2.l13=ds.campcen2.l13 + 1890
s.campcen.max.l13=ds.campcen.max.l13 + 1890
s.campcen.min.l13=ds.campcen.min.l13 + 1890


### Lecavalier et al., 2017: Camp Century Elevation change (12-0 kyr) ########
file=paste(fldr,"/lecavalier2017/datasets/Lecavalier-etal_Elevation_CampCent.tab",sep="")
data=read.table(file, skip=14, header=F)
time=-rev(data$V1) # kyr BP
time=time*1000 # yr BP
 
ds.campcen.l17=rev(data$V2)        # m 
ds.campcen.max.l17=rev(data$V3)    # 2 sigma
ds.campcen.min.l17=rev(data$V4)    # 2 sigma

s.campcen.l17=ds.campcen.l17 + 1890
s.campcen.max.l17=ds.campcen.max.l17 + 1890
s.campcen.min.l17=ds.campcen.min.l17 + 1890

t.ds.l17=time  # yr BP

### Lecavalier et al., 2017: Agassiz temperature anomaly (12 -0 kyr)  #####
file=paste(fldr,"/lecavalier2017/datasets/Lecavalier-etal_d18O_AgassizIceCap.tab",sep="")
data=read.table(file, skip=20, header=F)
time=-rev(data$V1)  # kyr BP
time=time*1000      # yr BP
tann=rev(data$V3)   # Mean annual T anom [°C]
tann.max=rev(data$V4) # Max annual anom T [°C]
tann.min=rev(data$V5) # Min annual anom T [°C]

tann.agas.l17=tann
tann.agas.max.l17 = tann.max
tann.agas.min.l17 = tann.min

t.tann.l17=time

### Lecavalier et al., 2017: Agassiz Summer temperature  (12-0 kyr) ####
file=paste(fldr,"/lecavalier2017/datasets/Lecavalier-etal_SummerT_AgassizIceCap.tab",sep="")
data=read.table(file, skip=18, header=F)
time=-rev(data$V1)  # kyr BP
time=time*1000      # yr BP
tsum=rev(data$V5)   # summer T [°C]
tsum.max=rev(data$V6) # summer T max [°C]
tsum.min=rev(data$V7) # min summer T [°C]

tsum.agas.l17=tsum
tsum.agas.max.l17=tsum.max
tsum.agas.min.l17=tsum.min

t.tsum.l17=time

### Kindler et al., 2014:  NGRIP T and accumulation #############
file=paste(fldr,"/kindler2014/kindler2014.csv",sep="")
data=read.csv(file, skip=12, sep=";", header=F)
time=-rev(data$V3)   # yr b2k
time=time[5656:length(time)]
time=time + 50 # yr BP
acc=rev(data$V4)
acc=acc[5656:length(acc)]  # m i.e./yr
tann=rev(data$V5)
tann=tann[5656:length(tann)]

acc.ngrip.k14=acc*0.917  # m i.e./yr -> m w.e./yr
tann.ngrip.k14=tann  # °C
t.k14=time # yr


##### Cuffey & Clow 1999: GISP2 accumulation #########
file=paste(fldr,"/grip-gisp/DATA/GISP2/PHYSICAL/ACCUM.DAT",sep="")
data=read.table(file, skip=51, sep="", nrows=10038)
time=-rev(data$V1)  # yr BP
acc=rev(data$V4) # m/yr i.e.

#ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
#acc=ma(acc, n=100)

acc.gisp2.cc99=acc*0.917 # m i.e. /yr -> m w.e./yr
t.cc99=time


##### Alley et al. 2000: GISP2 accumulation ##########
file=paste(fldr,"/alley2000/gisp2_accum_alley2000.txt",sep="")
data=read.table(file, skip=71, sep="", nrows=1768)
time=-rev(data$V1)*1000 # yr BP
acc=rev(data$V2)  # m ice /yr

acc.gisp2.a00=acc*0.917 # m i.e./yr -> m/yr w.e.
t.acc.a00=time



##### Alley et al. 2000: GISP2 annual temperature ######
file=paste(fldr,"/alley2000/gisp2_temp_alley2000.txt",sep="")
data=read.table(file, skip=75, sep="", nrows=1707)
time=-rev(data$V1)*1000 # yr BP
tann=rev(data$V2)  # °C

tann.gisp2.a00=tann
t.tann.a00=time


##### Kobashi et al., 2017: GISP2 Holocene temperature #####
file=paste(fldr,"/kobashi2017/kobashi2017-data.csv",sep="")
data=read.csv(file, skip=18, sep=";", header=F)
time=-rev(data$V10)  # yr BP
tann=rev(data$V11)   # °C
tann.max=rev(data$V12)  # °C
tann.min=rev(data$V13)  #°C

tann.gisp2.k17=tann
tann.max.gisp2.k17=tann.max
tann.min.gisp2.k17=tann.min
t.tann.k17=time


##### Rassmussen et al. 2013 (CP): NEEM accumulation #####
file=paste(fldr,"/rassmussen2013/acc-2013-12-05GICC05modelext-NEEM-1.csv",sep="")
data=read.csv(file, skip=2, sep=";", header=F)
time=-rev(data$V2) # yr b2k
time=time + 50 # yr BP
acc.mean=rev(data$V3) # m/yr i.e.
acc.min=rev(data$V4) # m/yr i.e. 
acc.max=rev(data$V5) # m/yr i.e.

acc.neem.mean.r13=acc.mean*0.917 # m/yr i.e. -> m/yr w.e.
acc.neem.min.r13=acc.min*0.917   # m/yr i.e. -> m/yr w.e.
acc.neem.max.r13=acc.max*0.917   # m/yr i.e. -> m/yr w.e.
t.r13=time


##### Badgeley et al., 2020 (CP): NGRIP, GRIP and Dye3 accumulation ########
file=paste(fldr,"/badgeley2020/not_used/Accumulation_Rate_Records_for_Moderate_Precipitation_Scenario_Raw.csv",sep="")
data=read.csv(file, skip=14, sep=",", header=F)
time.grip=-rev(data$V5) # yr BP
acc.grip=rev(data$V6) # kg/m2/yr
time.ngrip=-rev(data$V7)
acc.ngrip=rev(data$V8)
time.dye3=-rev(data$V9)
acc.dye3=rev(data$V10)

acc.grip=acc.grip*10^(-3)  # kg/m2/yr -> m/yr w.e. 
acc.ngrip=acc.ngrip*10^(-3) # kg/m2/yr -> m/yr w.e.
acc.dye3=acc.dye3*10^(-3) # kg/m2/yr -> m/yr w.e.

acc.grip.b20=acc.grip
t.grip.b20=time.grip
acc.ngrip.b20=acc.ngrip
t.ngrip.b20=time.ngrip
acc.dye3.b20=acc.dye3
t.dye3.b20=time.dye3



##### Gkinis et al., 2014: NGRIP accumulation ######
file=paste(fldr,"/gkinis2014/NorthGRIP_accu_DJmodel.tab",sep="")
data=read.table(file, skip=16, sep="")
time=-rev(data$V2) # kyr b2k
time=time*1000 +50 # yr BP
acc=rev(data$V3)  # m i.e./yr

acc.ngrip.g14=acc*0.917 # m i.e./yr -> m w.e./yr
t.g14=time


##### Hammer & Dahl Jensen 1999:  GRIP accumulation #######
#file="/home/titan/gwgk/gwgk005h/work/ngrip/scripts/data/grip-gisp/DATA/GRIP/PHYSICAL/GRIPACUM.DAT"
#data=read.table(file, skip=31, sep="")
#time=-rev(data$V4) # yr BP
#acc=rev(data$V2)  # g/cm2/yr
#acc=acc*10/917  # m i.e./yr
#
#grip.acc.hdj99=acc
#t.hdj99=time
#
#### -> not sure about the units


#### Buizert ez al., 2018: NEEM, NGRIP, GISP2, Dye3, Renland, Agassiz, Camp century, EGRIP annual and summer temperature  #####
file=paste(fldr,"/buizert2018/buizert2018-ice-core-sites.csv",sep="")
data=read.csv(file, skip=12, sep=";", header=F)
time=-rev(data$V1)  # yr BP
tann.neem=rev(data$V2) # °C
tann.ngrip=rev(data$V8)
tann.gisp2=rev(data$V14)
tann.dye3=rev(data$V20)
tann.ren=rev(data$V26)
tann.aga=rev(data$V32)
tann.campcen=rev(data$V44)
tann.egrip=rev(data$V50)

tsum.neem=rev(data$V5)
tsum.ngrip=rev(data$V11)
tsum.gisp2=rev(data$V17)
tsum.dye3=rev(data$V23)
tsum.ren=rev(data$V29)
tsum.aga=rev(data$V35)
tsum.campcen=rev(data$V47)
tsum.egrip=rev(data$V53)

tann.neem.b18=tann.neem
tann.ngrip.b18=tann.ngrip
tann.gisp2.b18=tann.gisp2
tann.dye3.b18=tann.dye3
tann.ren.b18=tann.ren
tann.aga.b18=tann.aga
tann.campcen.b18=tann.campcen
tann.egrip.b18=tann.egrip

tsum.neem.b18=tsum.neem
tsum.ngrip.b18=tsum.ngrip
tsum.gisp2.b18=tsum.gisp2
tsum.dye3.b18=tsum.dye3
tsum.ren.b18=tsum.ren
tsum.aga.b18=tsum.aga
tsum.campcen.b18=tsum.campcen
tsum.egrip.b18=tsum.egrip

t.b18=time





