\name{Fsiland.quantile}
\alias{Fsiland.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile computation for  spatial influence functions
}
\description{
The function computes the radius for which the integrated influence for the 
different landscape variables are equal to p.
}
\usage{
Fsiland.quantile(x, p = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class Fsiland resulting from function Fsiland()
}
  \item{p}{
a vector of probabilities
}
}

\value{
returned value is matrix. The number of rows corresponds to the number of landscape variables, and 
the number of columns to the number of probabilities in argument p.
}

\author{
O. Martin and F. Carpentier
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\donttest{
data(dataSiland)
data(landSiland)
resF=Fsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,sif="exponential")
Fsiland.quantile(x=resF,p=c(0.5,0.95))
}

}

