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
 $ L1:'data.frame':	12243 obs. of  2 variables:
  ..$ X: num [1:12243] 1075.98 1119.021 1162.06 1205.09 1248.138 ...
  ..$ Y: num [1:12243] 7988.14 7988.14 7988.14 7988.14 7988.14 ...
 $ L2:'data.frame':	2721 obs. of  2 variables:
  ..$ X: num [1:3615] 903.82 946.86 989.90 1032.94 3787.45 ...
  ..$ Y: num [1:3615] 7988.14 7988.14 7988.14 7988.14 7988.14 ...
}

\note{The landscape has been obtained by a simualtion. Names of the two landscape variables are L1 and L2.}

\examples{
data(landSiland)
names(landSiland)
#locations for the two landscape variables b1 and b2
plot(landSiland[[1]],col=2,pch=".")
points(landSiland[[2]],col=3,pch=".")
}
\keyword{datasets}
