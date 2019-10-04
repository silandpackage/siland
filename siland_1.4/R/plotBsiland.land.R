plotBsiland.land=function(x,land,data,var=1,landviz=F)
{
  if(class(x)!="Bsiland")
    stop("Argument x have to be an object of class Bsiland")
  if(length(var)>1)
    stop("length(var)>1: only one variable is allowed")
  if(class(land)[1] != "sf")
    stop("Argument land has to be an object of class sf")
  if(class(data) != "data.frame")
    stop("Argument data has to be an object of data.frame")
  #res=x
  border=x$border
  #requireNamespace("ggplot2")
  data=x$newdata
  distb=x$parambuffer[var]
  coef=x$coefficients
  
  #nb=str_sub(names(distb),6)
  nb=unlist(strsplit(names(distb),".",fixed=T))[2]
  cat(paste("Plot for landscape variable ",nb))
  cat("\n")
  print(distb)
  nbd=which(colnames(data)==nb)
  param=coef[which(names(coef)==nb)]
  data$Effect=data[,nbd]*param
  
  ggsf=st_as_sf(land)
  loc.sf=st_as_sf(data,coords = c("X","Y"))
  st_crs(loc.sf)<-st_crs(ggsf)$proj4string
  landsf=st_as_sf(land)
  
  mmin=min(data$Effect)
  mmax=max(data$Effect)
  X="X"
  Y="Y"
  xlim=range(data[,"X"])
  ylim=range(data[,"Y"])
  
  if (border==F)
  {
  buf=st_buffer(loc.sf ,distb)
  buf$Effect=data$Effect
  p=ggplot(data = loc.sf)+geom_sf(data=buf,aes_string(fill="Effect"),alpha=0.8,col="lightgrey")+ 
    geom_sf(size=0.8)+
    scale_fill_gradient2(low = "#0000FF", mid = "white", high = "#CC0000", midpoint = 0, limits = c(mmin, mmax)) + 
    theme_classic()+theme(panel.grid.major = element_line(colour = "white"), 
          panel.background = element_rect(fill = "white", colour = NA),axis.title = element_blank(),
          legend.title = element_blank(), legend.position = "bottom")
  
  #p=ggplot(data = loc.sf)+geom_sf(aes_string(fill="Effect"),alpha=0.8,col="lightgrey")+ geom_sf(size=0.8)
    
  
  }
  else
  { 
   parsf=st_intersects(loc.sf,landsf)
   parsf=landsf[unlist(parsf),]
   bufparsf=st_buffer(parsf,distb)
   bufparsf$Effect=data[,nbd]*param
   
   p=ggplot(bufparsf)+geom_sf(data=bufparsf,aes_string(fill="Effect"),alpha=0.8,col="lightgrey")+
     geom_sf(data=parsf,col="lightgrey",fill="white")+
     scale_fill_gradient2(low = "#0000FF", mid = "white", high = "#CC0000", midpoint = 0, limits = c(mmin, mmax)) + 
     theme(panel.grid.major = element_line(colour = "white"), 
           panel.background = element_rect(fill = "white", colour = NA),axis.title = element_blank(),
           legend.title = element_blank(), legend.position = "bottom")  
   
  }
  if (landviz==T)
    {
    nbl=which(colnames(landsf)==nb)
    landsf=landsf[landsf[[nbl]]==1,]
    if (border==T)
      { 
      p=ggplot(bufparsf)+ geom_sf(data=landsf)+geom_sf(data=bufparsf,aes_string(fill="Effect"),alpha=0.8,col="lightgrey")+geom_sf(data=parsf,col="lightgrey",fill="white")+
      scale_fill_gradient2(low = "#0000FF", mid = "white", high = "#CC0000", midpoint = 0, limits = c(mmin, mmax)) + 
      theme(panel.grid.major = element_line(colour = "white"), panel.background = element_rect(fill = "white", colour = NA),axis.title = element_blank(),legend.title = element_blank(), legend.position = "bottom")  
      }
    else
     { 
      p=ggplot(loc.sf)+geom_sf(data=landsf)+geom_sf(data=buf,aes_string(fill="Effect"),alpha=0.8,col="lightgrey")+geom_sf(size = 0.8)+ 
      scale_fill_gradient2(low = "#0000FF", mid = "white", high = "#CC0000", midpoint = 0, limits = c(mmin, mmax)) + 
      theme(panel.grid.major = element_line(colour = "white"), panel.background = element_rect(fill = "white", colour = NA),axis.title = element_blank(),legend.title = element_blank(), legend.position = "bottom")  
     }
    }
  
  return(p) 
}