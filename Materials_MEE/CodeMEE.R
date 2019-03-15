library(siland)

#Data loading with data() function
data(dataCmoth)
data(landCmoth)

#Data loading with data.gis() and land.gis() functions
dataCmoth=data.gis(dsn="./GIS",layer="dataCarpo",varname=c("Cmoth","trait"))
landCmoth1=land.gis(dsn="./GIS",layer="landCarpo",varname="OrgConv",landname = c("conv","org"),wd=40)
landCmoth2=land.gis(dsn="./GIS",layer="landCarpo",varname="Vine",landname = "vine",wd=40)
landCmoth=c(landCmoth1,landCmoth2)
names(landCmoth)

#Estimation with a gaussian SIF
res3G=siland(loc.model=Cmoth~trait,land=landCmoth,data=dataCmoth,sif="gaussian",test=T,family=gaussian)
summary(res3G)

res2G=siland(loc.model=Cmoth~trait,land=landCmoth[c(1,2)],data=dataCmoth,sif="gaussian",test=T,family=gaussian)
summary(res2G)

#Estimation on a simulated example
?landSiland
?dataSiland
data("dataSiland")
data("landSiland")
siland(loc.model= y~locvar+(1|Id), land=landSiland, data=dataSiland, sif="exponential", test=T, family="poisson")

res2=siland(loc.model= Cmoth~trait, land=landCmoth[c(1,2)],data=dataCmoth, sif="exponential", test=T,family=gaussian)
summary(res2)

siland.quantile(res2G,c(0.5,0.95))

plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=1)
plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=2)
plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=0)

plotsiland(res2G,landCmoth[c(1,2)],dataCmoth)
plotsiland.sif(res2G)

