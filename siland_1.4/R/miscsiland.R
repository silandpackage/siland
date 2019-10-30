
#Siland is a R package to estimate spatial infleunce of landscape variables.
#Copyright (C) 2017  Martin O. <olivier.martin@inra.fr>
#		     Carpentier F.	

# Siland is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


print.Fsiland<-function(x,...){
  cat("Model: ")
  print(as.formula(x$formula),showEnv=F)
  cat("\n")

  cat("Landscape variables: ")
  cat(colnames(x$landcontri))
  cat("\n \n")
  
  cat("Coefficients:\n")
  print(round(x$coefficients,3))
  if(x$modelType=="LMM" ||x$modelType=="GLMM")
    {
    cat("\nRandom effects:\n")
    if(x$modelType=="LMM"||x$modelType=="GLMM")
      print(x$rand.StdDev)
    }

  cat("\nstandard error: ")
  cat(x$sd.error)
  cat("\n")
  cat(paste("AIC: ",round(x$AIC,2), "\t", "AIC (no landscape): ",round(x$AIC0,2),sep=""))
  cat("\n")
  cat("(No landscape effect) p-value: ")
  if(x$pval0 < 1e-16) cat("<1e-16") else cat(x$pval0)
  
  
  #cat("\t(No landscape effect) p-value: ")
  #if(x$pval0==0) cat("<1e-16") else cat(x$pval0) 
  #cat("\n")
}




summary.Fsiland<-function(object,...)
{
  #data farame summary
  
  x<-object
  summaryx=summary(x$result)
  summaryx$call=x$formula
  
  cat("SIF parameters:\n")
  print(round(x$paramSIF,4))
  cat("\n")
  cat("-- Tests are given conditionnaly to the best SIF parameters --" )
  cat("\n")
  print(summaryx)
  
}


calcdist=function(data,land)
{
  #compute distance between observation location of data
  #and point source location of list land
  tmp=data[,c("X","Y")]
  Distobs=list()
  p=length(land)
  for(i in 1:p){
    Distobs[[i]]=t(apply(tmp,1,function(x){sqrt((x[1]-land[[i]][,"X"])^2+(x[2]-land[[i]][,"Y"])^2)}))
  }
  names(Distobs)=paste("Dist",1:p,sep="")
  return(Distobs)
}


calcscontri=function(distmoy,Distobs,w=1,sif="exponential")
{
  #compute contributions of sources that arise to each observation location.
  p=length(Distobs)
  n=nrow(Distobs[[1]])
  if (length(distmoy)==1){distmoy=rep(distmoy,p)}
  contri=list(NULL)
  if(sif=="exponential")
    for(i in 1:p) contri[[i]]= fdispE(Distobs[[i]],distmoy[i])
  
  if(sif=="gaussian")
    for(i in 1:p) contri[[i]]= fdispG(Distobs[[i]],distmoy[i])
  
  if(sif=="uniform")
    for(i in 1:p) contri[[i]]= fdispU(Distobs[[i]],distmoy[i])
 
  scontri=matrix(0,ncol=p,nrow=n)
  for(i in 1:p) scontri[,i]=apply(as.matrix(contri[[i]]),1,sum)
  
  scontric=scontri*w^2
  
  return(scontric)
}



FsilandMinusLoglik<-function(d,Dist,land,data,formula,sif,family)
{
  #options(warn=-1)
#print(d)
  #compute the minus loglikelihood for parameter  
  # of fis fucntion, that is the mean distance
  #data are local observations
  #land are list of landscape variables
  for(i in 1:length(d))
  {
    if(d[i]<0) 
      {
      mloglik=10^6
      return(mloglik)
    }
    if(d[i]>2000 ) 
    {
      mloglik=10^6
      return(mloglik)
    }
  }
  #w=min(dist(land[[1]][1:10,c("X","Y")]))
  w=sqrt(abs(diff(sort(unique(land[[1]][,1]))[1:2])*diff(sort(unique(land[[1]][,2]))[1:2])))
  
  #Dist=calcdist(data,land) 
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  #print(d)
  
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(glm(as.formula(formula),data=newdata,family=family), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
 # options(warn=0)
  invisible(return(mloglik))
}

FsilandMinusLoglikLMM<-function(d,Dist,land,data,formula,sif,family)
{
  options(warn=-1)
  #compute the minus loglikelihood for parameter  
  # of fis fucntion, that is the mean distance
  #data are local observations
  #land are list of landscape variables
  for(i in 1:length(d))
  {
    if(d[i]<0) 
    {
      mloglik=10^6
      return(mloglik)
    }
  }
  #w=min(dist(land[[1]][1:10,c("X","Y")]))
  w=sqrt(abs(diff(sort(unique(land[[1]][,1]))[1:2])*diff(sort(unique(land[[1]][,2]))[1:2])))
  
  Dist=calcdist(data,land) 
  #landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(lmer(as.formula(formula),data=newdata,REML=F), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
  
  options(warn=0)
  
  invisible(return(mloglik))
}

FsilandMinusLoglikGLMM<-function(d,Dist,land,data,formula,sif,family)
{
  options(warn=-1)
  #compute the minus loglikelihood for parameter  
  # of fis fucntion, that is the mean distance
  #data are local observations
  #land are list of landscape variables
  for(i in 1:length(d))
  {
    if(d[i]<0) 
    {
      mloglik=10^6
      return(mloglik)
    }
  }
  #w=min(dist(land[[1]][1:10,c("X","Y")]))
  w=sqrt(abs(diff(sort(unique(land[[1]][,1]))[1:2])*diff(sort(unique(land[[1]][,2]))[1:2])))
  
  Dist=calcdist(data,land) 
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  #landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif)
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(glmer(as.formula(formula),data=newdata,family=family), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
  options(warn=0)
  invisible(return(mloglik))
}






fdispE=function(x,dmoy){
  #compute density for exponential sif 
  alpha=dmoy/2
  (1/(2*pi*alpha^2))*exp(-x/alpha)   
}

fdispG=function(x,dmoy){
  #compute density for gaussian sif 
  alpha=dmoy/gamma(3/2)
  (1/(alpha^2*pi))*exp(-(x/alpha)^2)   
}

fdispU=function(x,dmoy){
  #compute density for uniform sif 
  alpha=3*dmoy/2
  w=x
  s=(x<=alpha)
  w[s]=1/(pi*alpha^2)
  w[!s]=0
  return(w)
  
}

fdispRE=function(r,dmoy){
  #compute radius density for exponential fis function
  alpha=dmoy/2
  (r/(alpha^2))*exp(-r/alpha)   
}

fdispRG=function(r,dmoy){
  #compute radius density for gaussian fis function
  alpha=dmoy/gamma(3/2)
  (2/alpha^2)*r*exp(-r^2/alpha^2)   
}

fdispRU=function(r,dmoy){
  #compute radius density for uniform fis function
  alpha=dmoy*3/2
  w=rep(0,length(r))
  s=r<=alpha
  w[s]=1/(pi*alpha^2)
  return(2*pi*r*w)
}


quantileE=function(q=0.9,dm,l=3000)
{  
  #Find quantile for radius distribution for exponential fis
  vv=seq(0,4*dm,length=l)
  pas=vv[2]-vv[1]
  cc=0
  for(i in 1:length(vv))
  {  
    tmp=fdispRE(vv[i],dm)*pas
    cc=cc+tmp
    if(cc>q)
    {
      resq=vv[i]
      return(resq)
    }
  }
}

quantileG=function(q=0.9,dm,l=3000)
{ 
  #Find quantile for radius distribution for gaussian fis
  vv=seq(0,4*dm,length=l)
  pas=vv[2]-vv[1]
  cc=0
  for(i in 1:length(vv))
  {  
    tmp=fdispRG(vv[i],dm)*pas
    cc=cc+tmp
    if(cc>q)
    {
      resq=vv[i]
      return(resq)
    }
  }
}

quantileU=function(dm,q=0.9,l=3000)
{  
  #Find quantile for radius distribution for uniform fis
  vv=seq(0,5*dm,length=l)
  pas=vv[2]-vv[1]
  cc=0
  for(i in 1:length(vv))
  {  
    tmp=fdispRU(vv[i],dm)*pas
    cc=cc+tmp
    #print(cc)
    if(cc>q)
    {
      resq=vv[i]
      return(resq)
    }
  }
}



  parens <- function(x) paste0("(",x,")")


leg.col <- function(colr, niv){
#add bar color scale for plotcontri  
  n <- length(colr)
  bx <- par("usr")
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  #box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
  #            bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 30)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = colr[i], border = colr[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(niv), max(niv)),
       yaxt = "n", ylab = "",
       xaxt= "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = 0.3)
  #axis(side = 4, las = 2, tick = FALSE, line = .25)
}


