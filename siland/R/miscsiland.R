
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


print.siland<-function(x,...){
  cat("Local model:\n")
  print(as.formula(x$loc.model),showEnv=F)
  cat("\n")
  cat("Landscape effect: ")
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
  cat("\t(No landscape effect) p-value: ")
  if(x$pval0==0) cat("<1e-16") else cat(x$pval0) 
}

#print<-function(x,...){
#UseMethod("print",x) }

#summary=function(object,...){
#UseMethod("summary") }




summary.siland<-function(object,...)
{
  #data farame summary
x<-object
  namesland=colnames(x$landcontri)
  nland=length(namesland)
  ncoef=length(x$coeff)
  if(!is.null(x$pval))
  {
    l=length(x$pval)
    #respval=round(x$pval,5)
    #respval[x$pval==0]="<1e-16"
    
    cat("Coefficients:\n")
    print(x$coefficients)
    
    if(x$modelType=="LMM" ||x$modelType=="GLMM")
    {
      cat("\nRandom effects:\n")
      if(x$modelType=="LMM"||x$modelType=="GLMM")
        print(x$rand.StdDev)
    }
    cat("\n")
    cat("pvalue (L.R. Test):\n")
    print(signif(x$pval,3))
    cat("\n")
    cat(paste("AIC: ",round(x$AIC,2), "\t", "AIC (no landscape): ",round(x$AIC0,2),sep=""))
    cat("\n")
    cat("(No landscape effect) p-value: ")
    if(x$pval0==0) cat("<1e-16") else cat(x$pval0)
  }
  else
  {
    cat("Coefficients:\n")
    print(x$coefficients)
    if(x$modelType=="LMM" ||x$modelType=="GLMM")
    {
    cat("\nRandom effects:\n")
    if(x$modelType=="LMM"||x$modelType=="GLMM")
      print(x$rand.StdDev)
    }
    cat("\n")
    cat(paste("AIC: ",round(x$AIC,2), "\t", "AIC (no landscape): ",round(x$AIC0,2),sep=""))
    cat("\n")
    cat("(No landscape effect) p-value: ")
    if(x$pval0==0) cat("<1e-16") else cat(x$pval0)
  }
}






#fonction caluculan entre le centre du verger et les lands qui ne se trouvent dans le verger
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


calcscontri=function(distmoy,Distobs,idland=NULL,idobs=NULL,w=1,sif="exponential")
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
  #contri<<-contri
  
  scontri=matrix(0,ncol=p,nrow=n)
  for(i in 1:p) scontri[,i]=apply(as.matrix(contri[[i]]),1,sum)
  
  # deletion of contributions from sources that are 
  # located in crops ds  les parcelles
  if(!is.null(idland) && !is.null(idobs))  
  {
    vv=vector("list",p)
    for(i in 1:p)
    {
      #for(j in 1:n)
      #{
      #  k=which(idland[[i]]%in%idobs[j])
      #  vv[[i]][j]=sum(as.matrix(contri[[i]])[j,k])
      #}
    }
    #for(i in 1:p)  
      #scontri[,i]=scontri[,i]-vv[[i]]
    scontric=scontri*w^2
  }
  else
  {
    scontric=scontri*w^2
  }
  return(scontric)
}



silandMinusLoglik<-function(d,data,land,formula,sif,family)
{
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
  w=min(dist(land[[1]][1:10,c("X","Y")]))
  Dist=calcdist(data,land) 
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  #landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif)
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(glm(as.formula(formula),data=newdata,family=family), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
  
  invisible(return(mloglik))
}

silandMinusLoglikLMM<-function(d,data,land,formula,sif,family)
{
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
  w=min(dist(land[[1]][1:10,c("X","Y")]))
  Dist=calcdist(data,land) 
  #landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(lmer(as.formula(formula),data=newdata,REML=F), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
  
  invisible(return(mloglik))
}

silandMinusLoglikGLMM<-function(d,data,land,formula,sif,family)
{
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
  w=min(dist(land[[1]][1:10,c("X","Y")]))
  Dist=calcdist(data,land) 
  landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif,w=w)
  #landcontri=calcscontri(distmoy=d,Distobs=Dist,sif=sif)
  colnames(landcontri)=names(land)
  newdata=cbind(data,landcontri)  
  if( inherits(rr <- try(glmer(as.formula(formula),data=newdata,family=family), silent = TRUE), "try-error"))
    mloglik= 10^6
  else  
    mloglik=as.numeric(-logLik(rr))
  
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


plotcontri<-function(res,land,data,type=0,numvar=NULL)
{
  #plot contributions of local variables and spatial variables
  #res: siland object
  #land: list of landscape associated to res
  #data : data file associated to res
  #type =type of plot, default type=0
  #numvar= numero of var for plot, by default numvar=NULL
  cexl=0.4
  cexo=0.8
  defpar=par()$mar
  par(mar=c(5.2,4,6,4))
  
  minx=min(unlist(lapply(land,function(x){min(x$X)})),data$X)
  miny=min(unlist(lapply(land,function(x){min(x$Y)})),data$Y)
  maxx=max(unlist(lapply(land,function(x){max(x$X)})),data$X)
  maxy=max(unlist(lapply(land,function(x){max(x$Y)})),data$Y)
  
  pp=c("X","Y")
  lr=length(land)
  nobs=nrow(data)
  if(class(res)!="siland")
    stop("argument res has to be of class siland")
  if(!is.numeric(type))
    stop("argument type has to be numeric")
  #if(type<0 |type>lr)
  #  stop("bad value for argument type")
  ### warning  if(is.null(land)) 
  
   #compute landscape contribution 
    landname=names(land) 
    coefland=matrix(rep(res$coefficients[landname],nobs),byrow=T,nrow=nobs)
    Contrilandscape=res$landcontri*coefland
    colnames(Contrilandscape)=landname
    
    #compute local contribution
    v=as.character(as.formula(res$loc.model)[[3]])
    localname=v[!(v=="+")]
    Contrilocal=NULL
    if(localname!="1")
    {
    coefloc=matrix(rep(res$coefficients[localname],nobs),byrow=T,nrow=nobs)
    Contrilocal=data[,localname]*coefloc
    colnames(Contrilocal)=localname
    }
    
    #All contributions
    contri=cbind(Contrilocal,Contrilandscape)
    Sumcontrilocal=NULL
    if(localname!="1") Sumcontrilocal=apply(Contrilocal,2,sum)
    Sumcontriland=apply(Contrilandscape,2,sum)
    Sumcontri=cbind(Sumcontrilocal,Sumcontriland)
    allname=landname
    if(localname!="1") allname=c(localname,landname)
    
    #Color graduation max and min
    #mmin=min(min(Sumlandscpe),min(Sumlocal))
    mmin=0
    mmax=max(abs(contri))  
    ncol=10
    
    #scale color
    col.gamme.pos=rainbow(ncol,start=4/6,end=1)
    col.gamme.neg=rainbow(ncol,start=2/6,end=4/6)
    col.gamme=c(col.gamme.neg,col.gamme.pos)
    scale.col.gamme=seq(-1,1,length.out=2*ncol)
    
    #red=0 yellow=1/6 green=2/6 cyan=3/6 blue=4/6 magenta=5/6
    contri=(contri-mmin)/(mmax-mmin)
        
    #number of plots
    nplot=ncol(contri)
    Mcolor=c()
    for (i in 1:nplot)
    { 
    color=rep(NA,nobs)
    posneg=contri[,i]<0
    color[posneg]= col.gamme.neg[ncol-ceiling(abs(contri[posneg,i])*(ncol-1))+1]
    pospos=contri[,i]>=0
    color[pospos]= col.gamme.pos[ceiling(contri[pospos,i]*(ncol-1))+1]
    Mcolor=cbind(Mcolor,color)
    }
    
  #type==0
  # contribution far all local and landscape variables are displayed
  if(type==0 & is.null(numvar))
  {
    #plot for local contribution
    if(!is.null(Contrilocal))
    {
    for (i in 1:ncol(Contrilocal))
      { 
      print(paste("local:",localname[i]))
      plot(data[,pp],pch=16,col=Mcolor[,i],main=paste(localname[i],": local contribution"),xlim=c(minx,maxx),ylim=c(miny,maxy),cex=cexo)
      leg.col(col.gamme, scale.col.gamme )
      readline(prompt = "Press Enter for the next graphic:\n")
    }
    }
    #plot for landscape contribution
    for (i in 1:ncol(Contrilandscape))
      {
      print(paste("spatial:",landname[i]))
      plot(land[[i]][,pp],col=grey(0.8),pch=16,xlim=c(minx,maxx),ylim=c(miny,maxy),main=paste(landname[i],": spatial contribution"),cex=cexl)
      if(!is.null(Contrilocal))
        points(data[,pp],pch=16,col=Mcolor[,i+ncol(Contrilocal)],cex=cexo)
      else
        points(data[,pp],pch=16,col=Mcolor[,i],cex=cexo)
      leg.col(col.gamme, scale.col.gamme)
      readline(prompt = "Press Enter for the next graphic:\n")
    }  
  }#end if type==0
  
  #plot for local contribution minus spatial contriibution
  if(type==1& is.null(numvar))
    {  
    if(localname=="1") stop("No local variable in the model")
    #diffcontri=as.vector(Sumcontri[,1]-Sumcontri[,2])
    diffcontri=as.vector(abs(Sumcontri[,2])/(abs(Sumcontri[,1])+abs(Sumcontri[,2])))
    diffcontri<<-diffcontri  
    mmax=max(abs(diffcontri))     
    ncol=10    
    #scale color
    col.gamme.pos=rainbow(ncol,start=4/6,end=1)
    col.gamme.neg=rainbow(ncol,start=2/6,end=4/6)
          
    #for a positive scale only
    col.gamme=c(col.gamme.pos)
    scale.col.gamme=seq(0,1,length.out=ncol)
    #scale with negative and positive values
    #col.gamme=c(col.gamme.neg,col.gamme.pos)
    #scale.col.gamme=seq(-mmax,mmax,length.out=2*ncol)
    ###red=0 yellow=1/6 green=2/6 cyan=3/6 blue=4/6 magenta=5/6
    diffcontri=(diffcontri)/(mmax)    
    color=rep(NA,nobs)
    posneg=diffcontri<0
    color[posneg]= col.gamme.neg[ceiling(abs(diffcontri[posneg])*(ncol-1))+1]
    pospos=diffcontri>=0
    color[pospos]= col.gamme.pos[ceiling(diffcontri[pospos]*(ncol-1))+1]    
    #print("Percent of spatial contribution")
    plot(data[,pp],type="n",col=color,main="Rate of landscape contribution",xlim=c(minx,maxx),ylim=c(miny,maxy))
    for(j in 1:length(land))
    {
      
    points(land[[j]][,pp],col=grey(0.8),pch=16,cex=cexl)
    }
    points(data[,pp],pch=16,col=color,cex=cexo)
    leg.col(col.gamme, scale.col.gamme )
    }
  
  #plot selectionned variables with numvar argument
  if(!is.null(numvar))
    {
      if(!is.vector(numvar) & !is.numeric(numvar))
        stop("varlocal argument has to be a numeric vector")
      if(max(numvar)> length(allname))
        stop("max of varlocal argument is too large")
      for (i in numvar)
      { 
        myvar=allname[i]
        #print(paste("local:",myvar))
        plot(data[,pp],pch=16,col=Mcolor[,i],main=paste(myvar," contribution"),xlim=c(minx,maxx),ylim=c(miny,maxy))
        leg.col(col.gamme, scale.col.gamme )
        if(length(numvar)>1){
        readline(prompt = "Please press Enter for the next graphic:\n")}
      }          
    }#end od if(!is.null(varlocal)) 
  par(defpar)
  invisible(return(NULL))
}

############################################################
#
#
#####################

plotsiland<-function(res,land,data)
{  
  q=0.9
  def.par <- par(no.readonly = TRUE)
  nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(0.8,0.2), c(0.2,0.8), TRUE)
  layout.show(nf)
  #landscape plot
  lx=range(unlist(lapply(land,function(x){range(x[,"X"])})))
  ly=range(unlist(lapply(land,function(x){range(x[,"Y"])})))
  
  pp=c("X","Y")
  par(mar = c(3,1,0,0))
  plot(land[[1]][,pp],xlim=lx,ylim=ly,type="n")
  for(i in 1:length(land ))
  {
    points(land[[i]][,pp],pch=16,cex=0.4,col=i+1 )
  }
  if(!is.null(data))
  {
    if(!is.null(data$y)){
      cmin=min(data$y)
      cmax=max(data$y)
      ccex=(data$y-cmin)/(cmax-cmin)
      points(data[,pp],pch=16,col=1,cex=ccex+0.2)
    }
    else
      points(data[,pp],pch=16,col=1,cex=1)
  }
  
  #plot of spatial inluence functions  
  sel=which(substring(names(res$coef),1,3)%in%"SIF")
  ldmoy=length(sel)
  dmoy=res$coef[sel]
  qsif=rep(0,ldmoy)
  for(i in 1:ldmoy)
  {
    if(res$sif=="exponential")
      qsif[i]=quantileE(q,dmoy[i])
    if(res$sif=="gaussian")
      qsif[i]=quantileG(q,dmoy[i])
    if(res$sif=="uniform")
      qsif[i]=quantileU(q,dmoy[i])
  }
  par(mar = c(0,1,2,0))
  plot(0,0,xlim=lx,ylim=c(-1,ldmoy+1),type="n",axe=F)
  axis(3)
  for(i in 1:ldmoy)
  {
    basex=lx[1]
    points(c(basex,basex+dmoy[i]),c(i,i),lwd=3,col=i+1,type="l")
    points(c(basex,basex+qsif[i]),c(i+0.3,i+0.3),lwd=3,lty=9,col=i+1,type="l")
    #print(dmoy[i])
  }
  
  # vertical plot
  par(mar = c(3,0,0,0))
  plot(0,0,ylim=ly,xlim=c(-1,ldmoy+1),type="n",axe=F)
  axis(4)
  for(i in 1:ldmoy)
  {
    basey=ly[1]
    #points(c(i,i),c(0,dmoy[i]),lwd=3,col=i+1,type="l")
    points(c(i,i),c(basey,basey+dmoy[i]),lwd=3,col=i+1,type="l")
    points(c(i+0.3,i+0.3),c(basey,basey+qsif[i]),lwd=3,lty=9,col=i+1,type="l")
    #print(dmoy[i])
  }
  par(def.par)
}



plotPays<-function(r,Dat=NULL,grid=NULL)
{
  if(is.null(Dat))
  {
  lx=range(unlist(lapply(r,function(x){range(x[,"X"])})))
  ly=range(unlist(lapply(r,function(x){range(x[,"Y"])})))
  }
  else
  {
    lxr=range(unlist(lapply(r,function(x){range(x[,"X"])})))
    lyr=range(unlist(lapply(r,function(x){range(x[,"Y"])})))
    lxd=range(Dat[,"X"])
    lyd=range(Dat[,"Y"])
    lx=range(c(lxr,lxd))
    ly=range(c(lyr,lyd))
  }
  
  pp=c("X","Y")
  if(is.null(grid))
    plot(r[[1]][pp],xlim=lx,ylim=ly,type="n",xlab="X",ylab="Y")
  else
    plot(0,0,xlim=grid,ylim=grid,type="n",xlab="X",ylab="Y")
  for(i in 1:length(r ))
  {
    points(r[[i]][,pp],col=i+1,pch="." )
  }
  if(!is.null(Dat))
  {
    if(!is.null(Dat$y)){
      cmin=min(Dat$y)
      cmax=max(Dat$y)
      ccex=(Dat$y-cmin)/(cmax-cmin)
      
      
      if(!is.null(Dat$locfac))
      {
        sel=Dat$locfac==1
        points(Dat[sel,pp],pch=16,col=1,cex=ccex[sel]+0.2)
        points(Dat[!sel,pp],pch=16,col=1,cex=ccex[!sel]+0.2)      
      }}
    else
      points(Dat[,pp],pch=16,col=1)
  }
}



rowsumsiland=function(x)
{
  if(is.vector(x)) {return(x)}else
  {return(apply(x,1,sum))}
}









