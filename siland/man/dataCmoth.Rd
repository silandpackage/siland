\name{dataCmoth}
\alias{dataCmoth}
\docType{data}
\title{
Data set for codling moth counts
}
\description{
The dataframe dataCmoth gives measurements for codling moth 
counts in 54 orchards. X and Y give the locations for observed data. The number Cmoth is the count of collected
 larvae in cardborad traps wrapped around tree trunks. trait is a local variable that gives the number of traitements for codling moth in the orchard.

}
\usage{data("dataCmoth")}
\format{
  A data frame with 54 observations on the following 4 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{Cmoth}}{a numeric vector}
    \item{\code{trait}}{a numeric vector}
    
  }
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
\keyword{datasets}
