landtoraster<-function(landgis,landname,wd,data=NULL)
{
  if(!is.null(data))
    border=T
  else
    border=F
  
  if(border==F)
  {
    lobs=list(NULL)
    extentLand=st_bbox(landgis)
    for (i in 1:length(landname))
    {
    
    r=raster(landgis,nrow=round((extentLand["ymax"]-extentLand["ymin"])/wd), ncol=round((extentLand["xmax"]-extentLand["xmin"])/wd),ext=extent(landgis))
    rland=fasterize(landgis,r,field=landname[i],fun='max')
    rlandpos=as.data.frame(rasterToPoints(rland))
    #only pixels different from zero are kept 
    rlandpos=rlandpos[rlandpos[,3]!=0,]
    colnames(rlandpos)=c("X","Y",landname[i])
    lobs[[i]]=rlandpos
    }
   resraster=lobs
  }
  
  if(border==T)
  {
    loc.sf=st_as_sf(as.data.frame(data),coords = c("X","Y"))
    st_crs(loc.sf)<-st_crs(landgis)$proj4string
    resraster=list(NULL)
    for(k in 1:length(landname))
    {
      tmpsig=landgis[landname[k]]
      stinter=st_intersects(loc.sf,landgis[landname[k]])
      geosel=unlist(lapply(stinter,function(x){if(length(x)==1) return(x) else return(-1000)} ))
      if(sum(geosel==-1000)>0)
        stop("Some observations are located outside of the boudaries of GIS landscape")
      tmpsig[geosel,][[landname[k]]]=0
      resraster[[k]]=landtoraster(tmpsig,landname = landname[k],wd=wd)[[1]]
    }
  }
  
  names(resraster)=landname
  return(landtable=resraster)
  
}
