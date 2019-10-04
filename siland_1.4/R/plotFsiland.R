plotFsiland=function(res,land,data)
{
  requireNamespace("ggplot2")
  requireNamespace("ggforce")
  namesland=colnames(res$landcontri)
  landraster=landtoraster(land,landname=namesland ,wd=res$wd)
  
  lx=range(unlist(lapply(landraster,function(x){range(x[,"X"])})))
  ly=range(unlist(lapply(landraster,function(x){range(x[,"Y"])})))
  
  pp=c("X","Y")
  Area=Fsiland.quantile(res,c(0.5,0.95))
  MR=max(Area)
  Area=rbind(0,Area)
  ldata=NULL
  DA=NULL
  
  
  xx=lx[2]+MR
  yy=ly[1]
  for(i in 1:length(landraster))
  {
    ldata=rbind(ldata,data.frame(landraster[[i]][,pp],LandVar=names(landraster)[i]))
    yy=yy+Area[i,2]+Area[i+1,2]+0.1*(ly[2]-ly[1])
    DA=rbind(DA,data.frame(x=xx,y=yy,r=Area[i+1,],LandVar=names(landraster)[i],Type=as.factor(c(0.5,0.95))))
  }
  
  cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  #version le rectangle sur les zones d'influence
 # g=ggplot(ldata,aes(x=X,y=Y,color=LandVar))+geom_point( )+
#    scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+
#    theme_classic()+theme(axis.title=element_blank(),legend.title=element_blank(),legend.position="bottom")+
#    geom_point(data=data,aes(x=X,y=Y),color=1)+coord_fixed()
#  g2=g+geom_rect(aes(xmin =lx[2],xmax= lx[2]+2*MR,ymin=ly[1],ymax=yy+Area[i+1,2]),fill="grey95",inherit.aes = FALSE)+
#    geom_circle(data=DA,aes(x0=x, y0=y, r=r),fill="White",colour="white",inherit.aes = FALSE,show.legend = FALSE)+
#    geom_circle(data=DA,aes(x0=x, y0=y, r=r,color=LandVar,fill=LandVar,linetype=Type,alpha=0.5),inherit.aes = FALSE,show.legend = FALSE)
#  g2
  
  #version le rectangle sur la zone d'?tude
  g=ggplot(ldata,aes_string(x="X",y="Y",color="LandVar"))+
    geom_rect(aes(xmin =lx[1],xmax= lx[2],ymin=ly[1],ymax=ly[2]),fill="grey95",inherit.aes = FALSE)+geom_point( size=0.4)+
    scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+
    theme_classic()+theme(axis.title=element_blank(),legend.title=element_blank(),legend.position="bottom")+
    geom_point(data=data,aes_string(x="X",y="Y"),color=1)+coord_fixed()
  g2=g+geom_circle(data=DA,aes_string(x0="x", y0="y", r="r",color="LandVar",fill="LandVar",linetype="Type",alpha=0.5),inherit.aes = FALSE,show.legend = FALSE)
  g2
  
  
}






