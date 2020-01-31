\name{bufferforsiland}
\alias{bufferforsiland}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes percentage of land use for buffers.
}
\description{
Computes percentage for landscape variables. This function is adapted when ladscape varibales correspond to land uses, such as conventional, orchard, forest, etc... 
}
\usage{
bufferforsiland(d, sfGIS, loc.sf, landnames, border = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
vector of distances to compute percentages for landscape uses given in landnames argument. The length of d has to be equal to that of landnames.
}
  \item{sfGIS}{
an object of class sf that gives the landscape variables (land uses). For land uses, landscape variables are equal to 1 or 0.
}
  \item{loc.sf}{
an object of class sf or a dataframe that gives locations for center buffers. In the case of a dataframe, locations have to be given in columns "X" and "Y".
}
  \item{landnames}{
vector of characters that gives the names of landscape variables for which buffer are computed
}
  \item{border}{
a logical (False by default). If False distance are computed from the point given in loc.sf. If True distance are computed from the boudnary of the polygon that contain the points in loc.sf.
}
}

\value{
a dataframe. Each column gives the percentage of landscape variables given in landnames and for buffer distances in d argument.
}
\author{
Martin O.
}


\examples{
\donttest{
data(dataSiland)
data(landSiland)
dataSiland[1:10,c("X","Y")]
landSiland$L1
landSiland$L2

bufferforsiland(c(100,200),sfGIS=landSiland, loc.sf=dataSiland, landnames=c("L1","L2"))

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
