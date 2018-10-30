\name{landCmoth}
\alias{landCmoth}
\docType{data}
\title{
Data sets for landscape description: organic and conventional features.
}
\description{
landCmoth is a list with three components. The first one gives the location of pixels for 
conventional variable, the second one for organic variable and the thrid for vine variable.
}
\usage{data("landCmoth")}
\format{
  The format is:
List of 3
$ conventional:'data.frame':	11282 obs. of  3 variables:
  ..$ X   : num [1:1099] 854982.5  855461  854373.5  854939  854982.5 ...
  ..$ Y   : num [1:1099] 6309549.4 6309549 6309512.6 6309513 6309512.6 ...
  ..$ Conv: num [1:1099] 1 1 1 1 1 1 1 1 1 1 ...
 $ organic     :'data.frame':	331 obs. of  3 variables:
  ..$ X  : num [1:36] 852372.4  852415.9  855852.5  855896  855939.5 ...
  ..$ Y  : num [1:36] 6307047.8 6307047.8 6307047.8 6307048 6307047.8 ...
  ..$ Bio: num [1:36] 1 1 1 1 1 1 1 1 1 1 ...
  $ vine      :'data.frame':	1977 obs. of  3 variables:
  ..$ X  : num [1:36] 851676.4  851763.4  851806.9  851850.4  851893.9 ...
  ..$ Y  : num [1:36]  6309365.5 6309365.5 6309365.5 6309365.5 6309365.5 ...
  ..$ Bio: num [1:36] 1 1 1 1 1 1 1 1 1 1 ...
 
}

\references{
Ricci B., Franck P., Bouvier J.-C., Casado D. and Lavigne C. (2011). Effects of hedgerow characteristics on intra-orchard distribution of larval codling moth. Agriculture, Ecosystems & Environment. 140. 395-400. 10.1016/j.agee.2011.01.001. 
}



\examples{
\donttest{
data(landCmoth)
data(dataCmoth)
#conventional pixels
plot(landCmoth[[1]][,c("X","Y")],pch=16,cex=0.3,col=1)
#organic pixels
points(landCmoth[[2]][,c("X","Y")],pch=16,cex=0.3,col=2 )
#vine pixels
points(landCmoth[[3]][,c("X","Y")],pch=16,cex=0.3,col=3 )
#measurement locations
points(dataCmoth[,c("X","Y")],pch=16,cex=0.8,col=4)
}
}

