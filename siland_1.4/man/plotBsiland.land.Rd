\name{plotBsiland.land}
\alias{plotBsiland.land}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial representation of the estimated buffers with Bsiland.
}
\description{
The function gives a spatial representation of buffers for the different variables. Representation can be done for each variable in the estimated model.
}
\usage{
plotBsiland.land(x, land, data, var = 1,landviz=F )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class Bsiland obtained with function Bsiland().
}
  \item{land}{
an object of class sf that gives the landscape variables for observations.
}
  \item{data}{
a dataframe containing the response variable and the local variables. The dataframe must be the same than the one used to obtain the object x.
}
  \item{var}{
a numeric value that gives the landscape variable for 
the spatial representation. 
}
  \item{landviz}{
If landviz=T (by default), polygons where variable var is observed are colored in gray.
}


}

\value{
an object of class ggplot
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Bsiland
}
\examples{
\donttest{
data(dataSiland)
data(landSiland)
res=Bsiland(obs~x1+L1+L2,land=landSiland,data=dataSiland,sif="exponential",family=gaussian)

plotBsiland.land(x=res,land=landSiland,data=dataSiland)
plotBsiland.land(x=res,land=landSiland,data=dataSiland,var=1,lw=20)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
