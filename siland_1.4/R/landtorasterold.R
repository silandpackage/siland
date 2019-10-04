landtorasterold<-function(landgis,landname,wd)
{
  #options(warn=-1)
  #if(!is.na(sum(as.numeric(landname)))) cat("Warning : landname have not to be a number")
  #options(warn=0)
  
  #landGIS=readOGR(dsn=dsn,layer=layer)
  #landGIS=as.data.frame(landgis[,landname])
  #Transformation to binary landscape variables 
  
  #m=c()
  
  #for(j in 1:length(landname))
  #{
  #  tmp=as.vector(landGIS[,varname])
  #  tmp[is.na(tmp)]=0
  #  vec=ifelse(tmp==landname[j],1,NA)
  #  if (sum(vec,na.rm=T)==0){print(paste("No category",landname[j],"in the",varname,"variable."))}
  #  m=cbind(m,vec)
  #}
  
  #landGIS@data=as.data.frame(m)
  #colnames(landGIS@data)=landname
  
  #discretisation of landscape variable 
  #creation of a list of matrix with discretized landscape variables
  lobs=list(NULL)
  extentLand=st_bbox(landgis)
  for (i in 1:length(landname))
  {
    #print(i)
    #if(is.null(extentLand))
    r=raster(landgis,ncol=round(extentLand["ymax"]-extentLand["ymin"])/wd, nrow=round(extentLand["xmax"]-extentLand["xmin"])/wd)
    #r=raster(landgis,res=50)
    # extent(r)=extent(as.vector(extentLand))
    #rland=rasterize(landgis,r,landname[i],fun='max')
    rland=fasterize(landgis,r,field=landname[i],fun='max')
    rlandpos=as.data.frame(rasterToPoints(rland))
    #only pixels different from zero are kept 
    rlandpos=rlandpos[rlandpos[,3]!=0,]
    colnames(rlandpos)=c("X","Y",landname[i])
    lobs[[i]]=rlandpos
  }
  
  names(lobs)=landname
  return(landtable=lobs)
  
}
