land.gis<-function(dsn,layer,varname,landname,wd=100,extentLand=NULL)
  #path.name  char describing path to folder containing data
  #layerland type char, name of gis layer file (shape file and database file)   without extension (.shp,.dbf)
  #containing variables landscape representing by polygons
  #varname type char, name of landscape variable presents in layerland.dbf
  #landname vector of names of labels of landscape variable. Labels have not to be numbers.
  #wd unit of discretization of polygons of landscape files
  #return a list of matrix of discretized surface location of each vallandvar 
  #extentLand : an object of class extent (see raster package)
{
  
  #data read
  uu=as.numeric(landname)
  if(!is.na(sum(as.numeric(landname)))) cat("Warning : landname have not to be a number")
  landSIG=readOGR(dsn=dsn,layer=layer)
  
  #Transformation to binary landscape variables 
  colnameslSIG=colnames(landSIG@data)
  m=c()
 
  for(j in 1:length(landname))
  {
    vec=c( )
    vec=ifelse(landSIG@data[,varname]==landname[j],1,NA)
    if (sum(vec,na.rm=T)==0){print(paste("No category",landname[j],"in the",varname,"variable."))}
    m=cbind(m,vec)
  }
  
  
  landSIG@data=cbind(landSIG@data,m)
  colnames(landSIG@data)=c(colnameslSIG,landname)
  
  #discretisation of landscape variable 
  #creation of a list of matrix with discretized landscape variables
  lobs=list(NULL)
  for (i in 1:length(landname))
  {
    if(is.null(extentLand))
      extentLand=extent(landSIG)
    r=raster(ncol=round(extent@ymax-extent@ymin)/wd, nrow=round(extent@xmax-extent@xmin)/wd)
    extent(r)=extentLand
    
    rland=rasterize(landSIG,r,landname[i],fun='max')
    rlandpos=as.data.frame(rasterToPoints(rland))
    colnames(rlandpos)=c("X","Y",landname[i])
    lobs[[i]]=rlandpos
  }
  names(lobs)=landname
  return(landtable=lobs)
}
