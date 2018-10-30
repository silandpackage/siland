siland.quantile<-function(x,p=0.95)
{

  if (class(x)!="siland"){
    stop("Argument x has to be of class siland")}
  if(!is.vector(p))
    stop("Argument p has to be a vector")
  if(sum(p>=1) |sum(p<0))
    stop("Argument p has to be a vector of probabilities (0<p<1) ")
  nl=ncol(x$landcontri)
  sif=x$sif
  namesLand=colnames(x$landcontri)
  coefsif=x$coefficients[namesLand]
  namesSIF=paste("SIF.",namesLand,sep="")
  valsif=x$coefficients[namesSIF]
  
  resq=matrix(0,ncol=length(p),nrow=length(valsif))
  for(i in 1:length(valsif))
  for(j in 1:length(p))
    {
      if(sif=="exponential")
        resq[i,j]=quantileE(p[j],valsif[i])
      if(sif=="gaussian")
        resq[i,j]=quantileG(p[j],valsif[i])
    }
    #resq=cbind(valsif,resq)
    colnames(resq)=as.character(p)
    rownames(resq)=namesSIF
    return(resq)
  }
  