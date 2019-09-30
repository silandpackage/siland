
#landCmoth=st_read(dsn = "~/Work/silandPackageGit/Materials_LandscapeEcology/GIS",layer="landCarpo")
#dataCmoth=st_read(dsn = "~/Work/silandPackageGit/Materials_LandscapeEcology/GIS/",layer="dataCarpo")

#org=ifelse(landCmoth$OrgConv=="org",1,0)
#org[is.na(org)]=0
#sum(org)
#conv=ifelse(landCmoth$OrgConv=="conv",1,0)
#conv[is.na(conv)]=0
#sum(conv)
#vine=ifelse(landCmoth$Vine=="vine",1,0)
#vine[is.na(vine)]=0
#landCmoth$conv=conv
#landCmoth$org=org
#landCmoth$vine=vine
#landCmoth=landCmoth[,-c(1,2)]

rm(list=ls())

library(siland)
library(sf)
library(lme4)
library(raster)
library(fasterize)
library(reshape2)
library(ggplot2)
library(ggforce)

landCmoth=st_read(dsn = "newCmoth",layer="landCmoth")
dataCmoth=st_read(dsn = "newCmoth",layer="dataCmoth")
dataCmoth=as.data.frame(cbind(X=dataCmoth$X,Y=dataCmoth$Y,trait=dataCmoth$trait,Cmoth=dataCmoth$Cmoth))


plot(landCmoth$geometry)
plot(landC[landCmoth$conv==1,]$geometry,col=2,add=T)
plot(landC[landCmoth$org==1,]$geometry,col=3,add=T)
plot(landC[landCmoth$vine==1,]$geometry,col=4,add=T)




res3G=Fsiland(Cmoth~trait+org+conv+vine,land=landC,data=dataC,sif="gaussian",init = c(50),wd=40)
res3G
summary(res3G)
plotFsiland.sif(res3G)
Fsiland.lik(resF,land=landC,data=dataC,varnames=c("org","conv","vine"),seqd=seq(5,2000,length=20))


res2G=Fsiland(Cmoth~trait+org+conv,land=landC,data=dataC,sif="gaussian",init = c(50),wd=40)
res2G

res2=Fsiland(Cmoth~trait+org+conv,land=landC,data=dataC,sif="exponential",init = c(100),wd=40)
res2


