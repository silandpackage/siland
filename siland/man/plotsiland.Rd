\name{plotsiland}
\alias{plotsiland}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot results from siland function
}
\description{
Plot locations of observed data and locations of landscape variables. The locations of black points correspond to the locations of observed data, and the size is propotional to numeric values. The locations of different landscape variables are displayed with small points of different colours. Estimated mean distance of spatial influence for each landscape variable is indicated beside the graphic with a continuous line. Dashed line indicates the radius that gives 95 percent of the total influence of each landscape variable. 
}
\usage{
plotsiland(res, land, data)
%%\method{print}{siland}(res,land,data)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{ an object obtained from the function siland
%%     ~~Describe \code{res} here~~
}
  \item{land}{a list with the location of the landscape variables. Each component of the list is a dataframe with two columns "X" and "Y" indicating the locations of each landscape variable.
%%     ~~Describe \code{land} here~~
}
  \item{data}{ a dataframe with the response variable and the local variables
%%     ~~Describe \code{data} here~~
}
}

\examples{
\donttest{
data(dataSiland)
data(landSiland)
res=siland(loc.model=y~locfac,land=landSiland,data=dataSiland,sif="exponential",family="gaussian")
plotsiland(res,landSiland, dataSiland)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
