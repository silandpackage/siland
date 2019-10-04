Fsiland<-function(formula,land,data,family="gaussian",sif="exponential",initSIF=50,border=F,wd=50)
{
  #land : a list of gis files (in case of several years) or one gis file (for one year) 
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
  
  
  #if(!is.list(land))
  #  stop("Problem: land argument has to be a list.")
  
  if(!is.data.frame(data))
    stop("data argument have to be a dataframe.")

  if(!(sif%in%c("gaussian","exponential","uniform")))
    stop("Problem: sif argument has to be \"gaussian\", \"exponential\" or \"uniform\" ")
  
  #if(is.null(names(land)))
   # stop("Problem: list for argument land has to have names ")
  
  
  model=as.formula(formula)
  modelType="GLM"
  termMix=findbars(model)
  if(!is.null(termMix))
  {
    if(isTRUE(all.equal(family,gaussian())))
      modelType="LMM"
    else
      modelType="GLMM"
  }
  
  if(class(land)[1]!="list")
    land=list(land)
  
  #sfGIS=lapply(land,st_as_sf)
  sfGIS=land

  if(class(data)[1]!="list")
    data=list(data)
  
  if(length(data)!=length(land))
    stop("The number of datasets for argument data have to be equal to the number of GIS objects for argument land")
  #test on column names for data
  datanames=colnames(data[[1]])
  if(length(data)>1)
  {
    for(i in 2:(length(data)))
    {
      tmpvari=colnames(data[[i]])
      if(sum(tmpvari==datanames)!=length(datanames))
        stop("Column names for data are not the same.")
    }
  }
  
  if(sum(datanames%in%c("X","Y"))!=2)
    stop("Error: colnames X and Y for observations are not available in data argument")
  
  landnames=names(sfGIS[[1]])

  allvars=all.vars(model)
  vary=allvars[1]
  
  
  #extract local variables
  allocalvars=datanames[datanames%in%allvars]
  localvars=allocalvars[!c(allocalvars%in%vary)]
  if(length(localvars)==0)
    localvars=NULL
  cat("Local variables: ")
  cat(paste(localvars[!c(localvars%in%vary)]," "))
  cat("\n")
  #extract landscape variables
  landvars=landnames[landnames%in%allvars]
  if(length(landvars)==0)
    stop("No landscape variable in the model")  
  cat("Landscape variables: ")
  cat(paste(landvars,sep=" "))
  cat("\n")
  if(length(localvars)==0)
    tmp=paste(landvars,collapse="+") 
  else
    tmp=paste(paste(localvars,collapse="+"), paste(landvars,collapse="+") , sep="+")
  model=as.formula(paste(vary, "~",tmp))
  
  cat("Model: ")
  print(model,showEnv=F)
  
  #model0 without landscape variables
  if(length(localvars[!c(localvars%in%vary)])>0)
    model0=as.formula( paste(paste(allvars[1],"~",sep=""),paste(localvars[!c(localvars%in%vary)],collapse="+"),sep="" ) )
  else
    model0=as.formula(paste(paste(allvars[1],"~",sep=""),paste("1",collapse="+"),sep="" ) )
  
  cat("Model0: ")
  print(model0,showEnv=F)
  
 
   # print(landvars)
  #initialisation algorithm
  if(is.null(initSIF))
    initSIF=rep(50,length(landvars))
  if(length(initSIF)==1)
    initSIF=rep(initSIF,length(landvars))
  if(length(initSIF)>1 & (length(initSIF)!=length(landvars)))
  stop("Problem : length(initSIF) and number of landscape variables have to be equal.")
  
  
  #transform dataframe data in sf object
  loc.sf=list(NULL)
  for(i in 1:length(data))
  {
    tmp=data[[i]][,c("X","Y")]
    loc.sf[[i]]=st_as_sf(tmp,coords = c("X","Y"))
    st_crs(loc.sf[[i]])<-st_crs(sfGIS[[i]])$proj4string
  }
  
 
 
  #compute raster from sig files.
  if(border==F)
  {
  landraster=landtoraster(sfGIS[[1]],landname = landvars,wd=wd)
  }

  if(border==T)
  {
  
    landraster=landtoraster(sfGIS[[1]],landname = landvars,wd=wd,data=data[[1]])
  }
    
  #landraster<<-landraster
  
  #surafce for a pixel in raster
  #suppose all pixels have same surface
  #w=min(dist(landraster[[1]][1:10,c("X","Y")]))
  w=sqrt(abs(diff(sort(unique(landraster[[1]][,1]))[1:2])*diff(sort(unique(landraster[[1]][,2]))[1:2])))
  
  DD=calcdist(data[[1]],landraster) 
  
  if(modelType=="GLM")
  {
    data1=data[[1]]
    myfun=function(d)
    {
    return(FsilandMinusLoglik(d,Dist=DD,data=data1,land=landraster,formula=model,sif=sif,family=family))
    }
  }
  if(modelType=="LMM")
  {
    data1=data[[1]]
    myfun=function(d)
    {
      return(FsilandMinusLoglikLMM(d,Dist=DD,data=data1,land=landraster,formula=model,sif=sif,family=family))
    }
  }
  if(modelType=="GLMM")
  {
    data1=data[[1]]
    myfun=function(d)
    {
      return(FsilandMinusLoglikGLMM(d,Dist=DD,data=data1,land=landraster,formula=model,sif=sif,family=family))
    }
  }
  
 
  if(length(landraster)>1)
  {
    #resoptim=optim(initSIF,myfun,lower=rep(0,length(land)),method="L-BFGS-B")
    #options(warn=-1)
    #ppar=rep(1/100,length=length(initSIF))
    #print(initSIF)
    resoptim=optim(initSIF,myfun)
    #resoptim=optim(par=initSIF,myfun,method="BFGS")
    #resoptim=nlm(myfun,initSIF)
    #options(warn=0)
  }
  
  if(length(landraster)==1)
  {
    resoptim=optimize(myfun,interval=c(0,1000))
    resoptim$par=resoptim$minimum
    options(warn=-1)
	#resoptim=optim(initSIF,myfun,method="Brent",lower=1,upper=2000)
	options(warn=0)
  }
  
  paramSIF=resoptim$par
  
  names(paramSIF)=paste("SIF.",landvars,sep="")
  
  landcontri=calcscontri(distmoy=paramSIF,Distobs=DD,sif=sif,w=w)
  colnames(landcontri)=landvars
  
  newdata=cbind(data[[1]],landcontri) 
  
  if(modelType=="GLM")
    restmp=glm(model,data=newdata,family=family ) 
  if(modelType=="LMM")
    restmp=lmer(model, data=newdata,REML=F ) 
  if(modelType=="GLMM")
    restmp=glmer(model,data=newdata,family=family ) 
  
  
  fit=predict(restmp)
  err=residuals(restmp)
  loglik=as.vector(logLik(restmp))
  sd.error=NA
  if(family$family=="gaussian" & modelType=="GLM")
    sd.error=sqrt(summary(restmp)$dispersion)
  if(family$family=="gaussian" &  modelType=="LMM")
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
    res0=glm(as.formula(model0),data=newdata,family=family)
    }
  if(modelType=="LMM")
    {
    res0=lmer(as.formula(model0),data=newdata,REML=F)
    }
  if(modelType=="GLMM")
    {
    res0=glmer(as.formula(model0),data=newdata,family=family)
    }
  loglik0=as.vector(logLik(res0))
  AIC0=as.vector(AIC(res0))
  
  # p.values for model with no landsacpe
  if(family$family!="gaussian")
    pval0=1-pchisq(2*(loglik-loglik0),2*length(landraster))
  else
    pval0=1-pchisq(2*(loglik-loglik0),2*length(landraster))
  if(modelType=="GLM"& family$family=="gaussian")
    nparam=length(resestim)+1
  if(modelType=="GLM" & family$family !="gaussian")
    nparam=length(resestim)
  if(modelType=="LMM" || modelType=="GLMM")
    nparam=length(resestim)+nrow(as.data.frame(rand.StdDev))
  
  
  resFsiland=list(coefficients=resestim,paramSIF=paramSIF, formula=model,landcontri=landcontri,loglik=loglik,loglik0=loglik0,result=restmp, fitted=fit,
                 sif=sif,resoptim=resoptim,result=restmp,AIC=AIC,AIC0=AIC0, nparam=nparam,pval0=pval0,family=family,
                 sd.error=sd.error,modelType=modelType,rand.StdDev=rand.StdDev,err=err,newdata=newdata,border=border,wd=wd)
  attr(resFsiland,"class") <- "Fsiland" 
  return(resFsiland)
}



