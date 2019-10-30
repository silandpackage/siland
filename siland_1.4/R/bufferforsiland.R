bufferforsiland<-function(d,sfGIS,loc.sf,landnames,border=F)
{
  #d is a vector 
  if(class(loc.sf)[1]!="sf")
  {
    loc.sf=st_as_sf(as.data.frame(loc.sf),coords = c("X","Y"))
    st_crs(loc.sf)<-st_crs(sfGIS)$proj4string
    
  }
  if(length(d) != length(landnames))
    stop("length(d) and length(landnames) have to be equal")
  
  n=nrow(loc.sf)
  nvar=length(landnames)
  #loc.sf=loc.sf[,c("X","Y")]
  #for each location in loc.sf, id object gives the polygon number in sfGIS
  #id=unlist(st_intersects(loc.sf,sfGIS))
    
  matBuffer=matrix(0,nrow=n,ncol=nvar)
  
  for(k in 1:nvar)
  {  
    
    #distrance are computed from observation locations
    if(!border) 
    {
    geom_buff=st_geometry(st_buffer(loc.sf,d[k]))
    #geom_buff0=st_geometry(st_buffer(sfGIS$geometry[loc2gis],0))
    #ll create a list with the geometry of geom_buff 
    ll=list(NULL)
    for(i in 1:length(geom_buff))
    {
      ll[[i]]=geom_buff[i]
    }
    
    }
    
    if(border)
    {
      stinter=st_intersects(loc.sf,sfGIS)
      loc2gis=unlist(lapply(stinter,function(x){if(length(x)==1) return(x) else return(-1000)} ))
      numerr=c(1:length(loc2gis))[loc2gis==-1000]
      if(sum(loc2gis==-1000)>0)
      {
        locerr=c(1:length(loc2gis))[loc2gis==-1000]
        cat("problem with observations : \n")
        print(locerr)
        stop("Some observations are not located inside polygons.")
      }
    geom_buff=st_geometry(st_buffer(sfGIS$geometry[loc2gis],d[k]))
    geom_buff0=st_geometry(st_buffer(sfGIS$geometry[loc2gis],0))
    #ll create a list with the geometry of geom_buff that does not intersect with geom_buff0
    ll=list(NULL)
    for(i in 1:length(geom_buff))
      ll[[i]]=st_difference(geom_buff[i],geom_buff0[i])
    }#end if(border)
    
    #areaBuff gives the area for the buffers in list ll
    areaBuff=unlist(lapply(ll,function(x){
                                          res=sum(st_area(x))
                                          if(length(res)==0)
                                            return(0)
                                          else
                                            return(res)
                                        }))
    
    #currebtland is the landscape for variable k 
    currentland=st_geometry(sfGIS[unlist(sfGIS[landnames[k]]) ==1,])
   
    #listArea gives the area for variable k and for the buffers in list ll
    listArea=lapply(ll,function(x){uu=st_intersection(x,currentland)
                                  res=as.numeric(sum(st_area(uu)))
                                  if(length(res)==0)
                                    return(0)
                                  else
                                    return(res)
                                  })  
    #transfomr listArea in vector
    Areavector=unlist(listArea)
    
    #propBuffer is proportion for variable k into all the buffers
    ind=c(areaBuff!=0)
    propBuffer=rep(0,length(areaBuff))
    propBuffer[ind]=Areavector[ind]/areaBuff[ind]
    matBuffer[,k]=propBuffer
  } #end for (k in 1:nvar)
 
   colnames(matBuffer)=landnames
   invisible(return(matBuffer)) 
}

