BsilandMinusLoglikLMM<-function(d,data,loc.sf,landnames,sfGIS,formula,family,border=F)
{
  options(warn=-1)
  
  for(i in 1:length(d))
  {
    if(d[i]<0) 
    {
      mloglik=10^6
      return(mloglik)
    }
  }
  
  
  matB=list(NULL)
  newdata=NULL
  for(i in 1:length(loc.sf))
  {
    matB[[i]]=bufferforsiland(d,sfGIS=sfGIS[[i]],loc.sf=loc.sf[[i]],landnames=landnames,border=border)
    newdata=rbind(newdata,cbind(data[[i]],matB[[i]]))
  }
  colnames(newdata)=c(colnames(data[[1]]),landnames)
  resout=lmer(formula,data=newdata)
  
  #if( inherits(rr <- try(glm(formula,data=newdata,family=family), silent = TRUE), "try-error"))
  #  mloglik= 10^6
  #else  
  mloglik=as.numeric(-logLik(resout))
  options(warn=0)
  invisible(return(mloglik))
}

  
  
  

