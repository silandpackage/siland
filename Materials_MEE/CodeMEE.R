library(siland)

#Data loading with data() function
data(dataCmoth)
data(landCmoth)


#Data loading with data.gis() and land.gis() functions
dataCmoth=data.gis(dsn="./GIS",layer="dataCarpo",varname="Cmoth")
landCmoth1=land.gis(dsn="./GIS",layer="landCarpo",varname="OrgConv",landname = c("Org","Conv"),wd=100)
landCmoth2=land.gis(dsn="./GIS",layer="landCarpo",varname="vine",landname = "vine",wd=100)
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
siland(loc.model= y~varloc+(1|Id), land=landSiland, data=dataSiland, sif="exponential", test=T, family="poisson")

res2=siland(loc.model= Cmoth~trait, land=landCmoth[c(1,2)],data=dataCmoth, sif="exponential", test=T,family=gaussian)
summary(res2)

siland.quantile(res2G,c(0.5,0.95))

plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=1)
plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=2)
plotsiland.land(res2G,land=landCmoth[c(1,2)],data=dataCmoth,var=0)



