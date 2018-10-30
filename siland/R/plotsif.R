plotsif<-function(d=NULL,sif="exponential" )
{
  #This function represents the density and the cumulative density 
  # for a given distance d and the type of spatial influence fucntion.
if(is.null(d)| (! is.numeric(d)))
  stop("argument d has to be a positive numeric value")
if(is.numeric(d) & d<0)
  stop("argument d has to be a positive numeric value")
  
if((sif != "exponential") & (sif!="gaussian"))
  stop("argument sif has to be \"gaussian\" or \"exponential\" ")
  
x=seq(0,4*d,length=200)
ua=x[2]
if(sif=="exponential")
{
  y=fdispE(x,d)
  yr=fdispRE(x,d)
}
if(sif=="gaussian")
{
  y=fdispG(x,d)
  yr=fdispRG(x,d)
}

par(mfrow=c(1,2))
plot(x,y,type="l",xlab="distance from source",ylab="density")
cyr=cumsum(yr)*ua
plot(x,cyr,type="l",xlab="radius around source",ylab="cumulative density")


return(invisible(list(density=cbind(x=x,y=y),cumdensity=cbind(x=x,y=cyr))))
}
