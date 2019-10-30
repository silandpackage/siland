\name{landSiland}
\alias{landSiland}
\docType{data}
\title{
A landscape with two landscape variables.
}
\description{
A landscape with two landscape variables, named L1 and L2.
}
\usage{data("landSiland")}
\format{
  A data frame with 4884 observations on the following 3 variables.
  \describe{
   
    \item{\code{geometry}}{a sfc_MULTIPOLYGON}
    \item{\code{L1}}{a numeric vector for landscape variable L1}
    \item{\code{L2}}{a numeric vector for landscape variable L2}
  }
}


\examples{
\donttest{
data(landSiland)
plot(landSiland$geometry)
plot(landSiland[landSiland$L1==1,]$geometry,col=2,add=T)
plot(landSiland[landSiland$L2==1,]$geometry,col=3,add=T)



}
}
\keyword{datasets}
