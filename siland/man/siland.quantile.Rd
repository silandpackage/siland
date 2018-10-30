\name{siland.quantile}
\alias{siland.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile computation for  spatial influence functions
}
\description{
The function computes the radius for which the integrated influence for the 
different landscape varaibles are equal to p.
}
\usage{
siland.quantile(x, p = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class siland resulting from function siland()
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
resE=siland(y~locvar,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)
siland.quantile(x=resE,p=c(0.5,0.95))
}

}

