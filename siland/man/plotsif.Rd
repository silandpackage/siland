\name{plotsif}
\alias{plotsif}

\title{
Plot density and cumultive density for sif function
}
\description{
The function represents the density and the cumulative density 
 for a given distance d and a shape for the spatial influence fucntion.

}
\usage{
plotsif(d = NULL, sif = "exponential")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
a positive numeric value that represents the parameter for the spatial influence function
}
  \item{sif}{
character that gives the shape of the spatial influence function, "exponential" (by default) or "gaussian".
}
}

\value{
a list of two dataframes. First element gives the density values for the spatial influence function, and the second gives the cumulative density.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
\donttest{
plotsif(d=300,sif="exponential")
}
}
