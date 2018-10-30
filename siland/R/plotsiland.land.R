plotsiland.land<-function(x,land, data, var=0,lw=100,xlim=NULL,ylim=NULL)
{
  #nsif=ncol(x$landcontri)
  #ll=length(x$coefficients)
  #valsif=x$coefficients[(ll-nsif+1):ll]
  w=min(dist(land[[1]][1:10,c("X","Y")]))
  rrx=unlist(lapply(land,function(x){range(x[,"X"])}))
  rry=unlist(lapply(land,function(x){range(x[,"Y"])}))
  
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
  nland=length(land)
  nl=length(x$coefficients)
  dres=x$coefficients[(nl-nland+1):(nl)]
  landsel=list(NULL)
  for(k in 1:length(land))
  {
    sel=(land[[k]][,"X"]>min(xr)) & (land[[k]][,"X"]<max(xr)) & (land[[k]][,"Y"]>min(yr)) & (land[[k]][,"Y"]<max(yr) )
    landsel[[k]]=land[[k]][sel,]
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
  mmax=max(c(abs(pcontri),abs(sumpcontri)))
  
  colorTable<- designer.colors(60, c( "darkblue","blue","white","red" ,"darkred") )
  #colorTable<- tim.colors(n=60,middle="white")
  brks<- seq(-mmax,mmax,length=61)
  
  
  if(var==0)
  {
    sumpcontri=apply(pcontri,1,sum)
    
    matcontri=matrix(sumpcontri,nrow=length(s.xr))
    image.plot(s.xr,s.yr,matcontri,breaks=brks,col=colorTable,xlab="X",ylab="Y")
    points(data[,c("X","Y")],pch=16,cex=0.8)
  }
    else
    {
     tmp=pcontri[,var] 
     matcontri=matrix(tmp,nrow=length(s.xr))
     image.plot(s.xr,s.yr,matcontri,breaks=brks,col=colorTable,xlab="X",ylab="Y")
     points(land[[var]][,c("X",("Y"))],pch=".",cex=1.2)
     points(data[,c("X","Y")],pch=16,cex=0.8)
    }
    
invisible(list(x=s.xr,y=s.yr,mat=matcontri))
}

