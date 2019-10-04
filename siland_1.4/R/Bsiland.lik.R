Bsiland.lik<-function(res,land,data,varnames=NULL,seqd=seq(2,2000,length=10))
{
  
  model=res$formula
  if(res$modelType=="GLMM")
    stop("This function don't take into account mixed effects")
  if(res$modelType=="LMM")
    stop("This function don't take into account mixed effects")
  
  seqd=sort(c(seqd,res$parambuffer))
  
  if(class(land)[1]!="list")
    land=list(land)
  sfGIS=lapply(land,st_as_sf)
  
  if(class(data)[1]!="list")
    data=list(data)
  
  if(length(data)!=length(land))
    stop("The number of datasets for argument data have to be equal to the number of GIS objects for argument land")
  
  loc.sf=list(NULL)
  for(i in 1:length(data))
  {
    tmp=data[[i]][,c("X","Y")]
    loc.sf[[i]]=st_as_sf(tmp,coords = c("X","Y"))
    st_crs(loc.sf[[i]])<-st_crs(sfGIS[[i]])$proj4string
  }
  
  datanames<-names(data[[1]])
  landnames=names(sfGIS[[1]])
  
  if(sum(varnames%in%landnames)<length(varnames))
    stop("Error: Some varnames are not variable in object land ")
  allvars=all.vars(model)
  vary=allvars[1]
  
  #extract landscape variables
  localvars=datanames[datanames%in%allvars]
  landvars=landnames[landnames%in%allvars]
  if(is.null(varnames))
    varnames=landvars
  
  nvars=length(varnames)
  matlik=matrix(0,nrow=nvars, ncol=length(seqd))
  selk= which(landvars %in% varnames)
  #print(varnames)
  for(k in 1:nvars)
  {
    cat(paste("Likelihood computing for ",varnames[k]))
    cat("\n")
    newd=res$parambuffer
    for(i in 1:length(seqd))
    {
      newd[selk[k]]=seqd[i]
      if(res$modelType=="GLM")
        matlik[k,i]= BsilandMinusLoglik(newd,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=res$family,border=res$border)
      if(res$modelType=="LMM")
        matlik[k,i]= BsilandMinusLoglikLMM(newd,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=res$family,border=res$border)
      if(res$modelType=="GLMM")
        matlik[k,i]= BsilandMinusLoglikGLMM(newd,data=data,loc.sf=loc.sf,landnames=landvars,sfGIS=sfGIS,formula=model,family=res$family,border=res$border)
    }
  }

  colnames(matlik)=seqd
  rownames(matlik)=varnames
  #matlik=cbind(id=1:nvars,t(matlik))
  matlik=rbind(matlik,"Estimated Model"=-res$loglik)
  matlik2=melt(matlik)
  colnames(matlik2)=c("Var1","Var2","value")
  linee=c("estimated model"=2)
  cols=c("estimated model"="darkorange")
  
 
  pp=ggplot(matlik2, aes_string(x="Var2",y="value",group="Var1"))+geom_line(aes_string(color="Var1"))+
    xlab("Buffer distance")+ylab("- Log-likelihood")+
    theme(legend.position="top",legend.title=element_blank())+
    geom_vline(xintercept=res$parambuffer[selk],color=1:nvars,lty=2)+
    geom_hline(yintercept=-res$loglik,lwd=1.5,color="darkorange")+
    scale_color_manual(values=c(1:(nvars),"darkorange"))
  #plot(pp)
    
  
  return(pp)
}


