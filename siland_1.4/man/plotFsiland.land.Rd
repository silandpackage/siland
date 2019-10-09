\name{plotFsiland.land}
\alias{plotFsiland.land}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial representation of the landscape influence for estimations with Fsiland() function
}
\description{
The function gives a spatial representation of landscape influence for the different variables. Representation can be done for each variable, or for the cumulative effects of all variables.
}
\usage{
plotFsiland.land(x, land, data, var = 0, lw = 100, xlim=NULL, ylim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class Fsiland obtained with function Fsiland().
}
  \item{land}{
an object of class sf that gives the landscape variables.
}
  \item{data}{
a dataframe containing the response variable and the local variables. The dataframe must be the same than the 
one used to obtain the object x.
}
  \item{var}{
a numeric value that gives the landscape variable for the spatial representation. If var=0 (by default), the  
representation gives the cumulative effects for all landscape variables.
}
  \item{lw}{
the number of pixels on x-axis and y-axis for the graphical representation. A too high value can lead to memory allocation problems. By default lw=100.
}
\item{xlim}{
a vector of two numbers (min and max on x-aixs) to give a representation for a sub-region. By default, xlim=NULL.
}
\item{ylim}{
a vector of two numbers (min and max on y-axis) to give a representation for a sub-region. By default, ylim=NULL.
}

}

\value{
an object of class ggplot
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Fsiland
}
\examples{
\donttest{
data(dataSiland)
data(landSiland)
res=Fsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)

plotFsiland.land(x=res,land=landSiland,data=dataSiland)
plotFsiland.land(x=res,land=landSiland,data=dataSiland,var=1,lw=20)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
