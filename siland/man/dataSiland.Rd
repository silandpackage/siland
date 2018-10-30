\name{dataSiland}
\alias{dataSiland}
\docType{data}
\title{
A simulated data set that reprsents observations.
}
\description{
The data set has 5 columns: X, Y, locfac, locvar and y. Columns X and Y represents the location for the observed data. Column locfac corresponds to an explanatory factor with two levels. Column locvar corresponds to a continuous explanatory variable. Column y corresponds to a conitnuous response variable simulated with a gaussian model. 
}
\usage{data("dataSiland")}
\format{
  A data frame with 80 simulated data and 5 variables.
  \describe{
    \item{\code{X}}{a numeric vector }
    \item{\code{Y}}{a numeric vector }
    \item{\code{locfac}}{a local factor with two levels}
    \item{\code{locvar}}{a continuous local variable}
    \item{\code{y}}{a numeric vector }
  }
}

\note{dataSiland$y have been simulated using the locfac variable and the landscape landSiland with a gaussian model. The locvar 
variable has not been used to simulate data y. To simulate the data, parameters have been fixed to the following values:
Intercept=-5, locfac=3, L1 =50, L2=-50, SIF.L1=400, SIF.L2=200 and sigma=2. The spatial influence function is exponential.

}

\examples{
data(dataSiland)
data(landSiland)
nrow(dataSiland)
#Plot for landscape variables
plot(landSiland[[1]],col=2,pch=".")
points(landSiland[[2]],col=3,pch=".")
#Locations of observations
points(dataSiland[,c("X","Y")])



}
\keyword{datasets}


