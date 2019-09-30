rm(list=ls())


library(siland)
library(sf)
library(lme4)
library(raster)
library(fasterize)
library(reshape2)
library(ggplot2)
library(ggforce)


data("dataSiland")
data("landSiland")

# Donnes d observation
dataSiland[1:10,]
# Donnees du paysage
landSiland[1:10,]
plot(landSiland$geometry)
plot(landSiland[landSiland$L1==1,]$geometry,col=2,add=T)
plot(landSiland[landSiland$L2==1,]$geometry,col=3,add=T)
points(dataSiland[,c("X","Y")],pch=16,col=4)


##############################"
##Estimation avec Fsiland
#############################
resF=Fsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,init = c(50),wd=20)
resF
summary(resF)

resF$landcontri

plotFsiland(resF,landSiland,dataSiland)

plotFsiland.land(resF,landSiland,dataSiland,var=2)
plotFsiland.land(resF,landSiland,dataSiland)

Fsiland.lik(resF,dataSiland,land=landSiland,varnames=c("L1","L2"),seqd=seq(50,500,length=20))

#Cas border=T (calucl des disntacesa partir du bord de la parcelle
resFT=Fsiland(obs~x1+L1+L2,data=dataSiland,land=landSiland,init = c(50),wd=20,border=T)
resFT
Fsiland.lik(resFT,dataSiland,land=landSiland,varnames=c("L1","L2"),seqd=seq(5,500,length=20))

#comparaison des contri enetre buffer=R et buffer=T
plot(resF$landcontri[,1],resFT$landcontri[,1])
abline(0,1)
plot(resF.F$landcontri[,2],resF.T$landcontri[,2])
abline(0,1)

#calcul de raster
landF=landtoraster(landSiland,c("L1","L2"),wd=20)
landFT=landtoraster(landSiland,c("L1","L2"),wd=20,data=dataSiland)

plot(landF[[1]][,1:2],pch=".")
points(landFT[[1]][,1:2],pch=".",col=2)
plot(landF[[2]][,1:2],pch=".")
points(landFT[[2]][,1:2],pch=".",col=2)

#cas de la sif unforme (un peu long)
resFu=Fsiland(obs~x1+L1+L2,data=dataSiland,land=landSiland,sif="uniform",init = c(50),wd=5)
resFu

##############################"
##Estimation avec Bsiland
#############################

resB=Bsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,initb = c(50))
resB
summary(resB)
plotBsiland.land(resB,landSiland,dataSiland,var=2)
resB$buffer

Bsiland.lik(resB,land=landSiland,data=dataSiland,varnames=c("L1","L2"),seqd=seq(50,1000,length=20))



bufferforsiland(c(100,200),sfGIS=landSiland, loc.sf=dataSiland, landnames=c("L1","L2"))


#Comparaison buffer et SIF uniforme
resB
resFu

plot(resB$buffer[,1],resFu$landcontri[,1] )
plot(resB$buffer[,2],resFu$landcontri[,2] )


####################################
#Cas de plusieurs annees avec Bsiland
####################################
##Calcul "un peu" long... 

year1=as.factor(rep("year1",nrow(dataSiland)))
year2=as.factor(rep("year2",nrow(dataSiland)))
dataSiland1=cbind("year"=year1,dataSiland)
dataSiland2=cbind("year"=year2,dataSiland)
dataSiland1[1:10,]
dataSiland2[1:10,]
ldataSiland=list(dataSiland1,dataSiland2)
llandSiland=list(landSiland,landSiland)

resB.y=Bsiland(obs~year+x1+L1+L2,land=llandSiland,data=ldataSiland,initb = c(50),border=F)

###

