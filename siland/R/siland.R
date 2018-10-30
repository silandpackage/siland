siland<-function(loc.model,land=NULL,data,initSIF=NULL,sif="exponential",family="gaussian",test=FALSE)
{
  #land : list of dataframe. Each dataframe gives the location ("X","Y") 
  #for the presence  of a landscape variable. The length of land list is the number of landscape varibale
  #data : a dataframe gives the variable of interest. Local variables have to be givent also in 
  #this datraframe. The location of the obervations have to be given by varible "X" and "Y"
  # check arguments
  

  if (is.character(family)) 
  {
    if(!(family%in%c("gaussian","poisson","binomial")))
      stop("Problem: family argument has to be \"gaussian\",\"poisson\" or \"binomial\" ")
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if(is.null(land))
  {
    stop("Problem: land argument has to be a non-null list.")
  }
  if(!is.list(land))
    stop("Problem: land argument has to be a list.")
  for(i in 1:length(land))
  {
    tmpnames=colnames(land[[i]])
    if(is.null(tmpnames) || is.na(tmpnames))
      stop(paste("Problem: No colnames in element",i,"of list land."))
    tmp=which(c("X","Y")%in%tmpnames)
    if(length(tmp)!=2)
      stop(paste("Problem: Missing X or Y colnames in element",i,"of list land."))
  }
  if(!is.data.frame(data))
    stop("data argument have to be a dataframe.")
  if(is.null(initSIF))
    initSIF=rep(100,length(land))
  if(length(initSIF)!=length(land))
    stop("Problem: length of initSIF argument is not equal to length of land.")
  if(!(sif%in%c("gaussian","exponential","uniform")))
    stop("Problem: sif argument has to be \"gaussian\", \"exponential\" or \"uniform\" ")
  
  if(is.null(names(land)))
    stop("Problem: list for argument land has to have names ")
  
  loc.model=as.formula(loc.model)
  modelType="GLM"
  termMix=findbars(loc.model)
  if(!is.null(termMix))
  {
    if(isTRUE(all.equal(family,gaussian())))
      modelType="LMM"
    else
      modelType="GLMM"
  }
  
  varmodel=all.vars(loc.model)
  vary=varmodel[1]
  namesland=names(land)
  #print(namesland)
  varx=NULL
  if(length(varmodel)>1)
  {
    varx=varmodel[2:length(varmodel)]
    if(sum(varx%in%colnames(data))!=length(varx))
      stop(paste("Problem: Error in loc.model specification. One explanatory variable is not available in data")) 
    if(sum(varx%in%namesland)>0)
      stop(paste("Problem: Local variables and landscape variables have to have different names.")) 
  }
  else
    varx=NULL
  
  if(!(vary%in%colnames(data)))
    stop(paste("Problem: Error in loc.model specification.",vary,"is not available in data"))

  #Minimum distance betwwen two points into the raster.
  w=min(dist(land[[1]][1:10,c("X","Y")]))
  
  landId=list(NULL)
  # for(i in 1:length(land))
  #  {
  #    landId[[i]]=NULL
  #    if(!is.null(land[[i]]$id))
  #      landId[[i]]=land[[i]]$id
  #  }
  
  formul=as.formula(paste(c(loc.model,namesland),collapse="+"))
  
  Dist=calcdist(data,land) 
  
  if(modelType=="GLM")
  {
    myfun=function(d)
    {
    return(silandMinusLoglik(d,data=data,land=land,formula=formul,sif=sif,family=family))
    }
  }
  if(modelType=="LMM")
  {
    myfun=function(d)
    {
      return(silandMinusLoglikLMM(d,data=data,land=land,formula=formul,sif=sif,family=family))
    }
  }
  if(modelType=="GLMM")
  {
    myfun=function(d)
    {
      return(silandMinusLoglikGLMM(d,data=data,land=land,formula=formul,sif=sif,family=family))
    }
  }
  
  if(length(land)>1)
  {
    #resoptim=optim(initSIF,myfun,lower=rep(0,length(land)),method="L-BFGS-B")
    resoptim=optim(initSIF,myfun)
  }
  
  if(length(land)==1)
  {
    #resoptim=optimize(myfun,interval=c(0,300000))
    #resoptim$par=resoptim$minimum
	resoptim=optim(initSIF,myfun)
  }
  
  paramSIF=resoptim$par
  names(paramSIF)=paste("SIF.",namesland,sep="")
  landcontri=calcscontri(distmoy=paramSIF,Dist,sif=sif,idland=NULL,idobs=NULL,w=w)
  colnames(landcontri)=namesland
  
  newdata=cbind(data,landcontri)   
  if(modelType=="GLM")
    restmp=glm(formul,data=newdata,family=family ) 
  if(modelType=="LMM")
    restmp=lmer(formul, data=newdata,REML=F ) 
  if(modelType=="GLMM")
    restmp=glmer(formul,data=newdata,family=family ) 
  
  fit=predict(restmp)
  err=residuals(fit)
  loglik=as.vector(logLik(restmp))
  sd.error=NA
  if(family$family=="gaussian" & modelType=="GLM")
    sd.error=sqrt(summary(restmp)$dispersion)
  if(family$family=="gaussian" & modelType=="LMM")
    sd.error=summary(restmp)$sigma
  rand.StdDev=NULL
  if(modelType=="GLM")
  {
  coeflm=restmp$coef
  resestim=c(coeflm,paramSIF)
  }
  if(modelType=="LMM"||modelType=="GLMM")
  {
    coeflm=summary(restmp)$coef[,1]
    resestim=c(coeflm,paramSIF)
    rand.StdDev=VarCorr(restmp)
  }
  
  nparam=length(resestim)
  if(family$family!="gaussian")
    AIC=2*nparam-2*loglik
  else
    AIC=2*(nparam+1)-2*loglik
  
  #Model with only local variables
  if(modelType=="GLM")
    {
    res0=glm(as.formula(loc.model),data=newdata,family=family)
    }
  if(modelType=="LMM")
    {
    res0=lmer(as.formula(loc.model),data=newdata,REML=F)
    }
  if(modelType=="GLMM")
    {
    res0=glmer(as.formula(loc.model),data=newdata,family=family)
    }
  loglik0=as.vector(logLik(res0))
  AIC0=as.vector(AIC(res0))
  
  # p.values for model with no landsacpe
  if(family$family!="gaussian")
    pval0=1-pchisq(2*(loglik-loglik0),2*length(land))
  else
    pval0=1-pchisq(2*(loglik-loglik0),2*length(land))
  if(modelType=="GLM"& family$family=="gaussian")
    nparam=length(resestim)+1
  if(modelType=="GLM" & family$family !="gaussian")
    nparam=length(resestim)
  if(modelType=="LMM" || modelType=="GLMM")
    nparam=length(resestim)+nrow(as.data.frame(rand.StdDev))
  
    #as.data.frame(rr$rand.StdDev)
  #compute pvalue for parameters by likelihood ratio test
  pval.land=rep(NULL,length(land))
  #pval.local=rep(NULL,length(varx))
  pval=NULL
  if (test==TRUE)
  {
    cat("-- Pvalues computing -- \n")
    #test for landscape variables
    if(length(land)==1)
    {
      pval.land=1-pchisq(2*(loglik-loglik0),2)
      names(pval.land)=namesland
    }
    
    if(length(land)>1)
    {
      for (i in 1:length(land))
      {
        res.test=siland(loc.model = loc.model,land=land[-i],data=data,sif=sif,initSIF=initSIF[-i],family=family,test=F)
        pval.land[i]=1-pchisq(2*(loglik-res.test$loglik),2)
      }
      names(pval.land)=namesland
    }
    
    #test for local variables
   
    #Bintercept="(Intercept)"%in%names(resestim)
    tmp=as.character(loc.model)[3]
    Lterms=unlist(strsplit(tmp,"\\+"))
    if(length(Lterms)==1)
    {
      modeltmp=as.formula(paste(vary, "1", sep="~" ))
      restmp=siland(loc.model = modeltmp,land=land,data=data,sif=sif,initSIF=initSIF,family=family,test=F)
      pval.local=1-pchisq(2*(loglik-restmp$loglik),nparam-restmp$nparam)
    }
  
    if(length(Lterms)>1)
    {
      pval.local=rep(0,length(Lterms))
      for(k in 1:length(Lterms))
      {
      modeltmp=as.formula(paste(vary, paste(Lterms[-k],collapse="+"), sep="~" ))
      restmp=siland(loc.model = modeltmp,land=land,data=data,sif=sif,initSIF=initSIF,family=family,test=F)
      pval.local[k]=1-pchisq(2*(loglik-restmp$loglik),nparam-restmp$nparam)
      #print(modeltmp)
      }
    }
   
    names(pval.local)=Lterms
    
    pval=c(pval.local,pval.land)
    
  }#end compute pvalues for parameters (if test ==TRUE)
  
  ressiland=list(coefficients=resestim,loc.model=loc.model,landcontri=landcontri,loglik=loglik,loglik0=loglik0,fitted=fit,
                 sif=sif,resoptim=resoptim,AIC=AIC,AIC0=AIC0, nparam=nparam,pval0=pval0,pval=pval,family=family,
                 sd.error=sd.error,modelType=modelType,rand.StdDev=rand.StdDev,nparam=nparam,err=err)
  attr(ressiland,"class") <- "siland" 
  return(ressiland)
}



