

plotFsiland.sif=function(x)
{
  requireNamespace("ggplot2")
  #if (class(x) != "siland") {
  #  stop("Argument x has to be of class siland")
  #}
  nl = ncol(x$landcontri)
  sif = x$sif
  namesLand = colnames(x$landcontri)
  coefsif = x$coefficients[namesLand]
  namesSIF = paste("SIF.", namesLand, sep = "")
  valsif = x$coefficients[namesSIF]
  fismax = max(valsif)
  xl = seq(0, 1. * fismax, length = 200)
  x2=seq(0, 1.8 * fismax, length = 1000)
  ua = xl[2]
  par(mfrow = c(nl, 2))
  density.fis = NULL
  bdensity.fis = NULL
 # data=data.frame(d=x2)
  data=NULL
  for (i in 1:nl) {
    if (sif == "exponential") {
      y = fdispE(x2, valsif[i])
      density.fis = cbind(density.fis, y)
      bdensity.fis = cbind(bdensity.fis, coefsif[i] * y)
      data=rbind(data,data.frame(d=x2,y=y,z= coefsif[i] * y,Landvar=namesLand[i]))
     
      
    }
    if (sif == "gaussian") {
      y = fdispG(x2, valsif[i])
      density.fis = cbind(density.fis, y)
      bdensity.fis = cbind(bdensity.fis, coefsif[i] * y)
      #data=cbind(data,y)
      #colnames(data)[i+1]=namesLand[i]
      data=rbind(data,data.frame(d=x2,y=y,z= coefsif[i] * y,Landvar=namesLand[i]))
      }
  }


  
  cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ggplot(data,aes_string(x="d",y="y",col="Landvar",fill="Landvar"))+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+
    geom_area(alpha=0.3,size=0.6)+theme_classic()+
    geom_vline(data=data.frame(dsif=valsif,Landvar=namesLand), aes_string(xintercept="dsif", color="Landvar"),size=1,show.legend = F)+
    labs(y="Spatial Influence Function",x="Distance",group="Landscape variable") +
    theme(legend.title=element_blank(),legend.justification=c(1,1), legend.position=c(1,1))
   
  # ggplot(data,aes_string(x=d,y=z,col=Landvar,fill=Landvar))+
  #   scale_x_continuous(expand = c(0, 0)) + 
  #   scale_fill_manual(values=cbPalette)+scale_color_manual(values=cbPalette)+
  #   geom_area(alpha=0.2,size=0.6)+theme_classic()+
  #   geom_vline(data=data.frame(dsif=valsif,Landvar=namesLand), aes_string(xintercept=dsif, color=Landvar),size=1,show.legend = F)+
  #   labs(y="Density",x="Distance",group="Landscape variable") +
  #   theme(legend.title=element_blank(),legend.justification=c(1,1), legend.position=c(1,1))
  
}



