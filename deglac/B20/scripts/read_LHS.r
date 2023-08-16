# LHS ensemble

fldr="/home/titan/gwgk/gwgk005h/work/ngrip/v1.801/8km/lhs/B18/deglac/B20"

# read file txt
filename = paste(fldr,"/lhs_np15_ns3000_values.txt", sep="")
txt=read.table(filename,header=F, sep="", skip=1)
a=as.matrix(txt)
colnames(a)=c("itmb","itmc","q","z0","nffdlt","enh","k","tau","lith","kt","negis_c","negis_s","negis_n","fp","negis1")    # order of columns depending on Yelmo output
nr=nrow(a)
nc=ncol(a)

itmb=round(as.vector(a[,1]),digits=4)
itmc=round(as.vector(a[,2]),digits=4)
q = round(as.vector(a[,3]),digits=4)
z0 = round(as.vector(a[,4]),digits=4)
neff = round(as.vector(a[,5]),digits=4)
enh = round(as.vector(a[,6]),digits=4)
k= round(as.vector(a[,7]),digits=4)
tau=round(as.vector(a[,8]),digits=4)
lith=round(as.vector(a[,9]),digits=4)
kt=round(as.vector(a[,10]),digits=4)
negis_c = round(as.vector(a[,11]),digits=4)
negis_s = round(as.vector(a[,12]),digits=4)
negis_n = round(as.vector(a[,13]),digits=4)
fp = round(as.vector(a[,14]),digits=4)
negis1=round(as.vector(a[,15]),digits=4)

itmb=sprintf('%.4f',itmb)
itmc=sprintf('%.4f',itmc)
q=sprintf('%.4f',q)
z0=sprintf('%.4f',z0)
neff=sprintf('%.4f',neff)
enh=sprintf('%.4f',enh)
k=sprintf('%.4f',k)
tau=sprintf('%.4f',tau)
lith=sprintf('%.4f',lith)
kt=sprintf('%.4f',kt)
negis_c=sprintf('%.4f',negis_c)
negis_s=sprintf('%.4f',negis_s)
negis_n=sprintf('%.4f',negis_n)
fp=sprintf('%.4f',fp)
negis1=sprintf('%.4f',negis1)

# find the simulation which has the parameters closest to the mean param values of the whole set of variables
m=c()
for(i in 1:nc){
  m[i]=median(a[,i])
}

min=rep(1e5,nc)
for(i in 1:nr){
 j=abs(a[i,]-m)
  if(sqrt(sum(j^2)) < sqrt(sum(min^2))){
    min[]=j[]
    median.par=i
  }
}


#min=rep(1e5,nc)
#for(i in 1:nr){
#  j=abs(a[i,]-m)
#  if(j[1]<min[1] & j[2]<min[2] & j[3]<min[3] & j[4]<min[4] & j[5]<min[5]){
#   min[]=j[]
#   median.par=i
#  }
#}

median.sim=a[median.par,]

# pass a, nr, nc, mean.sim


