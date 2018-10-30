\name{landSiland}
\alias{landSiland}
\docType{data}
\title{
A list of simulated data sets that describes landscape.
}
\description{
Each component of the list is a dataframe that correponds to a landscape variable. For each data set, the columns X and Y indicate the locations of the different "pixels".
}
\usage{data("landSiland")}
\format{
  The format is:
List of 2
 $ L1:'data.frame':	9513 obs. of  2 variables:
  ..$ X: num [1:9513] 0 43 86.1 129.1 172.2 ...
  ..$ Y: num [1:9513] 9535 9535 9535 9535 9535 ...
 $ L2:'data.frame':	3615 obs. of  2 variables:
  ..$ X: num [1:3615] 387 430 473 516 560 ...
  ..$ Y: num [1:3615] 9535 9535 9535 9535 9535 ...  
}

\note{The landscape comes from the north of France. Two features of the landscape are considered and they are denoted L1 and L2.}

\examples{
data(landSiland)
names(landSiland)
#locations for the two landscape variables b1 and b2
plot(landSiland[[1]],col=2,pch=".")
points(landSiland[[2]],col=3,pch=".")
}
\keyword{datasets}
