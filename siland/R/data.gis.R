data.gis<-function(dsn,layer,varname,factor.var=NULL)
{
  #dsn  char describing path to folder containing data
  #layer type char, name of gis layer files (shape file and database file)  without extension (.shp,.dbf)
  #varname type char, vector of names of  variables in layerdata.dbf file
  #factor.var vector of boolean specify wether variables has to be considered as factors or not
  #return a data.frame with column names x,y and varname v
  
  obsGIS=readOGR(dsn=dsn,layer=layer) 
  if(! (is.logical(factor.var) |is.null(factor.var)) ) 
    stop("Argument as.factor.localvar has to be NULL or logical")
  
  print(varname)
  if(sum(is.element(colnames(obsGIS@data),varname))!=length(varname))
    stop("Argument varname: some variables are not in OGR data source")
  
  dobs=as.data.frame(cbind(obsGIS@coords,obsGIS@data[,varname]))
  
  for (i in (1:length(varname))[factor.var])
  {
    dobs[,2+i]=as.factor(dobs[,2+i])
    #print(as.factor(dobs[,3+i]))
  }
    colnames(dobs)=c("X","Y",varname)  
  
  
  return(dobs)
}
