plotFsiland.land<-function(x,land, data, var=0,lw=100,xlim=NULL,ylim=NULL)
{
  requireNamespace("ggplot2")
  namesland=colnames(x$landcontri)
landraster=landtoraster(land,landname=namesland ,wd=x$wd)
  
  #for(i in 1:length(landraster))
  #{
  #  landraster[[i]]<-as.data.frame(landraster[[i]])
  #}
  w=min(dist(landraster[[1]][1:10,c("X","Y")]))
  rrx=unlist(lapply(landraster,function(x){range(x[,"X"])}))
  rry=unlist(lapply(landraster,function(x){range(x[,"Y"])}))
  
  xr=range(rrx)
  yr=range(rry)
  
  if(!is.null(xlim))
    xr=xlim
  if(!is.null(ylim))
    yr=ylim
  
  lx=xr[2]-xr[1]
  ly=yr[2]-yr[1]
  wd=round(max(lx/lw,ly/lw),0)
  
  s.xr=seq(xr[1],xr[2],by=wd)
  s.yr=seq(yr[1],yr[2],by=wd)
  
  mgrid=expand.grid(s.xr,s.yr)
  colnames(mgrid)<-c("X","Y")
  nland=length(landraster)
  nl=length(x$coefficients)
  dres=x$coefficients[(nl-nland+1):(nl)]
  landsel=list(NULL)
  for(k in 1:length(landraster))
  {
    sel=(landraster[[k]][,"X"]>min(xr)) & (landraster[[k]][,"X"]<max(xr)) & (landraster[[k]][,"Y"]>min(yr)) & (landraster[[k]][,"Y"]<max(yr) )
    landsel[[k]]=landraster[[k]][sel,]
  }
  print("Distance computing... Wait...")
  Dist=calcdist(mgrid,landsel)
  print("Contribution computing... Wait..")
  landcontri=calcscontri(distmoy=dres,Distobs=Dist,sif=x$sif,w=w)
  #landocntri is multiplied by the strength of the 
  #different landscape variables
  coef=x$coefficients[(nl-2*nland+1):(nl-nland)]
  nr=nrow(landcontri)
  mcoef= matrix(rep(coef,nr),nr,byrow=T)
  pcontri=mcoef*landcontri
  #pcontri=landcontri
  sumpcontri=apply(pcontri,1,sum)
  #mmax=max(c(abs(pcontri),abs(sumpcontri)))
  mmax=max(c(pcontri,sumpcontri))
  mmin=min(c(pcontri,sumpcontri))
  
  #colorTable<- designer.colors(60, c( "darkblue","blue","white","red" ,"darkred") )
  #colorTable<- tim.colors(n=60,middle="white")
 # brks<- seq(-mmax,mmax,length=61)
  
  if(var==0)
  {
    sumpcontri=apply(pcontri,1,sum)
    
    dd=cbind(sumpcontri,mgrid)
    g= ggplot(data=dd, aes_string(x="X", y="Y", fill = "sumpcontri")) + geom_raster(interpolate = T)+
      scale_fill_gradient2(low="#0000FF",mid="white",high="#CC0000",midpoint=0,limits=c(mmin,mmax))+
      coord_fixed()+
      theme_classic() + theme(axis.title=element_blank(),legend.title=element_blank(),legend.position="bottom") + 
      geom_point(data=data,aes_string(x="X",y="Y"),inherit.aes=F)
    
    }
    else
    {
      dd=cbind(V=pcontri[,var],mgrid)
      g= ggplot(data=dd, aes_string(x="X", y="Y"))+ geom_raster(interpolate = T,aes_string(fill ="V"))+
        scale_fill_gradient2(low="#0000FF",mid="white",high="#CC0000",midpoint=0,limits=c(mmin,mmax))+
        coord_fixed()+
        theme_classic()+theme(axis.title=element_blank(),legend.title=element_blank(),legend.position="bottom")+
        geom_point(data=landraster[[var]],aes_string(x="X",y="Y"),colour="#999999",fill="#999999",size=0.3,inherit.aes=F)+
        geom_point(data=data,aes_string(x="X",y="Y"),inherit.aes=F)
        
       
      
    }
    
  return(g)
}

