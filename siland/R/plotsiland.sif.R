plotsiland.sif<-function(x)
{

  #x must be an object of class siland
  #this function shows the density for the 
  #different estimated influence functions. 
  #The graphics give information on the form of the decrease 
  #for the different influence functions. These are not the cumulative 
  #effects of the landscape, but the effect of one pixel of the raster.
  
  if (class(x)!="siland"){
    stop("Argument x has to be of class siland")}
  nl=ncol(x$landcontri)
  sif=x$sif
  namesLand=colnames(x$landcontri)
  coefsif=x$coefficients[namesLand]
  namesSIF=paste("SIF.",namesLand,sep="")
  valsif=x$coefficients[namesSIF]
  fismax=max(valsif)
  
  xl=seq(0,1.8*fismax,length=200)
  ua=xl[2]
  par(mfrow=c(nl,2))
  
  ##densoty computing
  density.fis=NULL
  bdensity.fis=NULL
  
  for(i in 1:nl)
  {
    if(sif=="exponential")
    {
      y=fdispE(xl,valsif[i])
      density.fis=cbind(density.fis,y)
      bdensity.fis=cbind(bdensity.fis,coefsif[i]*y)
    }
    if(sif=="gaussian")
    {
      y=fdispG(xl,valsif[i])
      density.fis=cbind(density.fis,y)
      bdensity.fis=cbind(bdensity.fis,coefsif[i]*y)
    }
  }
  
  par(mfrow=c(1,2))
  matplot(xl,density.fis,type="l",xlab="distance from source",ylab="density",lty=1,col=1:nl)
  legend("top",legend=namesLand,col=1:nl,cex=0.8,lty=1,lwd=2)
  matplot(xl,bdensity.fis,type="l",xlab="distance from source",ylab="landscape coefficient * density",lty=1,col=1:nl)
  legend("top",legend=namesLand,col=1:nl,cex=0.8,lty=1,lwd=2)
  abline(h=0)
  
}