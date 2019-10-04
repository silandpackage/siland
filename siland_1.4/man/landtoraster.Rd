\name{landtoraster}
\alias{landtoraster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform an object class sf into raster points
}
\description{
The function transforms an object of class sf into raster points with wd the pixel sizes of raster. Object 
landmane gives the variable names for which raster points have to be computed.
}
\usage{
landtoraster(landgis, landname, wd, data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{landgis}{
an object of class sf 
}
  \item{landname}{
  Variable names for which raster points have to be computed.
}
  \item{wd}{
a number that correpond to size  pixels for raster
}
  \item{data}{
a dataframe with positions X and Y. If a dataset is specified, the pixels in polygons 
containing the points in the dataframe are deleted.
}
  
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a list of dataframes. The components of the list correspond to variables specifed in argument landname. Each dataframe gives the position X and Y for pixels in raster.
}


\author{
O. Martin
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\donttest{
data(dataSiland)
data(landSiland)
dataSiland[1:10,c("X","Y")]
landSiland$L1
landSiland$L2
u1=landtoraster(landSiland,c("L1","L2"),wd=20)
names(u1)
plot(u1[[1]][,c("X","Y")],pch=".")
u2=landtoraster(landSiland,c("L1","L2"),wd=20,data=dataSiland)
points(u2[[1]][,c("X","Y")],pch=".",col=2)
points(dataSiland[,c("X","Y")],pch=16,cex=0.6,col=3)


}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
