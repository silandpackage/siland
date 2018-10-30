\name{plotcontri}
\alias{plotcontri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot contributions
}
\description{
Plot contributions for the different landscape variables
}
\usage{
plotcontri(res, land, data, type = 0, numvar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
a result form siland estimation
}
  \item{land}{
the landscape associated to the estimation res
}
  \item{data}{
the observations associated to the estimation res
}
  \item{type}{
type of plot, 0 or 1. For type equal to 0 (by default), a graphic for each local and landscape variable is displayed. For type equal to 1, the rate of landscape with local contribtions is displayed.
}
  \item{numvar}{
  a number indicating the variable that has to be displayed. If NULL (by default), graphics are displayed following the argument type.
}
}
\note{
A local model with only continuous variables can be currently considered. The function doesn't work if local model includes random effect(s) or categorical variables (factor). These cases are ongoing  developements.  
}

\examples{
\donttest{
data(landSiland)
data(dataSiland)
res=siland(y~locvar,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)
plotcontri(res,landSiland,dataSiland)

plotcontri(res,landSiland,dataSiland,type=1)
plotcontri(res,landSiland,dataSiland,numvar=2)
}
}


