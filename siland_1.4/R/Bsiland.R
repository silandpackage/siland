Bsiland<-function(formula,land,data,family="gaussian",init=200,border=F)
{
  #formula : an object of class formula
  #land : a list of gis files (in case of several years) or one gis file (for one year)
  #for the presence  of a landscape variable. The length of land list is the number of landscape varibale
  #data : a data frame for observations
  
  # check arguments

  if (is.character(family)) 
  {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
 
  
  if(class(land)[1]!="list")
    land=list(land)
  sfGIS=lapply(land,st_as_sf)
  
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
  
  
  ##### extract local and landscapes variables
  landnames=names(sfGIS[[1]])
  termMix=findbars(model)
  

  
  #labels(terms(model))
  
  allvars=all.vars(model)
  allvars=all.vars(model)
  vary=allvars[1]
  allvars=allvars[-1]
  
  #localvars=datanames[datanames%in%allvars]
  #landvars=landnames[landnames%in%allvars]
  
  localvars=allvars[allvars%in%datanames]
  landvars=allvars[allvars%in%landnames]
  
  cat("Local variables: ")
  cat(paste(localvars,sep=" "))
  cat("\n")
  #extract landscape variables
  if(length(landvars)==0)
    stop("No landscape variable in the model")  
  cat("Landscape variables: ")
  cat(paste(landvars,sep=" "))
  cat("\n")
  
  
  if(length(localvars)==0)
    model0=as.formula(paste(vary,"~1"))
  labelsFormula=labels(terms(nobars(model)))
  
  sel=NULL
  for(j in 1:length(labelsFormula))
  {
    finteraction=sum(grep(":",labelsFormula[j]))
    if(finteraction==0)
      sel=c(sel,labelsFormula[j]%in%datanames)  
    else
    {
      u=unlist(strsplit(labelsFormula[j],":"))
      if(sum(u%in%landnames)>0)
        sel=c(sel,F)
      else
        sel=c(sel,T)
    }
  }
  
  if(is.null(termMix))
  {
    if(length(localvars)>0)
      model0=as.formula( paste(paste(vary,"~",sep=""),paste(labelsFormula[sel],collapse="+") ,sep="" ) )
    else
      model0=as.formula(paste(paste(vary,"~",sep=""),paste("1",collapse="+"),sep="" ) )
  }
  
  if(!is.null(termMix))
  {
    if(length(localvars)>0)
      model0=as.formula( paste(paste(vary,"~",sep="") ,paste(labelsFormula[sel],collapse="+"),"+",paste(parens(termMix),collapse="+")  ,sep=" " ) )
    else
      model0=as.formula( paste(paste(vary,"~",sep="") ,paste(parens(termMix),collapse="+")  ,sep=" " ) )
  }
  
  cat("Model: ")
  print(model,showEnv=F)
  
  cat("Model0: ")
  print(model0,showEnv=F)
  
  
  ##### end extrcat local and landscape variables 
  
  #initialisation algorithm
  if(is.null(init))
    init=rep(50,length(landvars))
  if(length(init)==1)
    init=rep(init,length(landvars))
  if(length(init)>1 & (length(init)!=length(landvars)))
stop("Probelem : length(init) and length(landnames) has to be equal.")
  
  
  
  #transform dataframe data in sf object
  loc.sf=list(NULL)
  for(i in 1:length(data))
  {
  tmp=data[[i]][,c("X","Y")]
  loc.sf[[i]]=st_as_sf(tmp,coords = c("X","Y"))
  st_crs(loc.sf[[i]])<-st_crs(sfGIS[[i]])$proj4string
  }
  
  
  if(modelType=="GLM")
  {
    myfun=function(d)
    {
    return(BsilandMinusLoglik(d,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=family,border=border))
    }
  }
  
  if(modelType=="LMM")
  {
    myfun=function(d)
    {
      return(BsilandMinusLoglikLMM(d,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=family,border=border))
    }
  }
  
  if(modelType=="GLMM")
  {
    myfun=function(d)
    {
      return(BsilandMinusLoglikGLMM(d,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=family,border=border))
    }
  }
  
  if(length(landvars)>1)
  {
    #resoptim=optim(initSIF,myfun,lower=rep(0,length(land)),method="L-BFGS-B")
    resoptim=optim(init,myfun)
  }
  
  if(length(landvars)==1)
  {
    #resoptim=optimize(myfun,interval=c(0,300000),minimum=T)
    #resoptim$par=resoptim$minimum
	resoptim=optim(init,myfun,method="Brent",lower=1,upper=5000)
  }
  
  paramBuffer=resoptim$par
 
  names(paramBuffer)=paste("B.",landvars,sep="")
  
  matB=list(NULL)
  newdata=NULL
  for(i in 1:length(loc.sf))
  {
    matB[[i]]=bufferforsiland(paramBuffer,sfGIS=sfGIS[[1]],loc.sf=loc.sf[[i]],landnames=landvars,border=border)
    newdata=rbind(newdata,cbind(data[[i]],matB[[i]]))
  }
  valBuffer=newdata[,landvars]
  if(modelType=="GLM")
    restmp=glm(model,data=newdata,family=family)
  if(modelType=="LMM")
    restmp=lmer(model,data=newdata)
  if(modelType=="GLMM")
    restmp=glmer(model,data=newdata,family=family)  
  
  
  fit=predict(restmp)
  err=residuals(restmp)
  loglik=as.vector(logLik(restmp))
  sd.error=NA
  if(family$family=="gaussian" & modelType=="GLM")
    sd.error=sqrt(summary(restmp)$dispersion)
  if(family$family=="gaussian" & modelType=="GLMM")
    sd.error=summary(restmp)$sigma
  rand.StdDev=NULL
  if(modelType=="GLM")
  {
    resestim=c(restmp$coef ,paramBuffer)
  }
  
  if(modelType=="LMM"||modelType=="GLMM")
  {
    resestim=c(summary(restmp)$coef[,1] ,paramBuffer)
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
    pval0=1-pchisq(2*(loglik-loglik0),2*length(landvars))
  else
    pval0=1-pchisq(2*(loglik-loglik0),2*length(landvars))
  if(modelType=="GLM" & family$family=="gaussian")
    nparam=length(resestim)+1
  if(modelType=="GLM" & family$family !="gaussian")
    nparam=length(resestim)
  if(modelType=="LMM" || modelType=="GLMM")
    nparam=length(resestim)+nrow(as.data.frame(rand.StdDev))
  
  #pval.local=rep(NULL,length(varx))
  #pval=NULL
  #coefestim=c(coefficients(res),paramBuffer)
  
  resBsiland=list(coefficients=resestim ,parambuffer=paramBuffer, formula=model,buffer=valBuffer,loglik=loglik,loglik0=loglik0,fitted=fit,
                  resoptim=resoptim,result=restmp,AIC=AIC,AIC0=AIC0, nparam=nparam,pval0=pval0,family=family,
                  sd.error=sd.error,modelType=modelType,rand.StdDev=rand.StdDev, err=err,newdata=newdata, border=border)
  attr(resBsiland,"class") <- "Bsiland" 
  return(resBsiland)
}


  
  