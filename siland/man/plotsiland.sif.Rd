\name{plotsiland.sif}
\alias{plotsiland.sif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the estimated spatial influence functions.
}
\description{
Graphic representations for the different influence functions estimated with the function siland().
}
\usage{
plotsiland.sif(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 must be an object of type siland, that is an object resulting from function siland().
}
}

\value{
no value is returned
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
siland
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\donttest{
data(dataSiland)
data(landSiland)
resE=siland(y~locvar,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)
plotsiland.sif(x=resE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
