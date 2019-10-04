\name{dataSiland}
\alias{dataSiland}
\docType{data}
\title{
A dataframe for simulated data
}
\description{
A dataframe for simulated observations and local variables.
}

\usage{data("dataSiland")}
\format{
  A data frame with 98 observations on the following 5 variables.
  \describe{
    \item{\code{X}}{location on X-axis}
    \item{\code{Y}}{location on Y-axis}
    \item{\code{x1}}{a continous local variable}
    \item{\code{Id}}{a numeric vector}
    \item{\code{obs}}{the observed variable}
  }
}


\examples{
\donttest{
data(dataSiland)
plot(dataSiland[,c("X","Y")])
}

}
\keyword{datasets}
