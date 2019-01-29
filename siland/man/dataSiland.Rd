\name{dataSiland}
\alias{dataSiland}
\docType{data}
\title{
A simulated data set that reprsents observations.
}
\description{
The dataframe has 5 columns: X, Y, y, locvar and Id. Columns X and Y represents the location for the observed data. Column y corresponds to a response variable simulated with a poisson model. Column locvar corresponds to a continuous explanatory variable. Id variable gives the number for the 40 statistcal units and each unit have 3 repeated measurements. The dataframe has 120 rows.
}
\usage{data("dataSiland")}
\format{
  A data frame with 120 simulated data and 5 variables.
  \describe{
    \item{\code{X}}{a numeric vector }
    \item{\code{Y}}{a numeric vector }
    \item{\code{y}}{a numeric vector }
    \item{\code{locvar}}{a continuous local variable}
\item{\code{Id}}{a numeric vector}
    
  }
}

\note{dataSiland$y have been simulated using the locvar variable and the landscape landSiland with a poisson model. Id variable gives the number for the 40 statistcal units and each unit have 3 repeated measurements. To simulate the data, parameters have been fixed to the following values:
Intercept=-1, locvar=1, L1 =-2, L2=-2, SIF.L1=100, SIF.L2=50. A random gaussian effect  (mean equal to zero and standard deviation equal to 0.5) has been added for observiations arising for a given statistical unit (Generalized Linear Mixed Model with poisson distrubtion). The spatial influence function is exponential.

}

\examples{
\donttest{
data(dataSiland)
data(landSiland)
nrow(dataSiland)
#Plot for landscape variables
plot(landSiland[[1]],col=2,pch=".")
points(landSiland[[2]],col=3,pch=".")
#Locations of observations
points(dataSiland[,c("X","Y")])
res=siland(y~locvar,land=landSiland,data=dataSiland,family=poisson)
summary(res)
#Model with random effect
res=siland(y~locvar+(1|Id),land=landSiland,data=dataSiland,family=poisson)
summary(res)
}
}
\keyword{datasets}


