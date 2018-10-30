\name{data.gis}
\alias{data.gis}

\title{
Read OGR vector maps for observed data or measurements.
}
\description{
The function reads an OGR data source and return a dataframe. The function is developped to 
read easily observed data or ecological measurements from OGR data source. A similar function (land.gis) 
is available to read OGR data source that describes the landscape associated to oberved data.
}

\usage{
data.gis(dsn, layer, varname, factor.var = NULL)
}

\arguments{
  \item{dsn}{ data source name (interpretation varies by driver, for some drivers, 
dsn is a file name, but may also be a folder)
}
  \item{layer}{layer name (varies by driver, may be a file name without extension). From rgdal 1.2.*, layer may be missing, in which case ogrListLayers examines the dsn, and fails if there are no layers, silently reads the only layer if only one layer is found, and reads the first layer if multiple layers are present, issuing a warning that layer should be given explicitly.
}
  \item{varname}{
a vector with the variable names that need to be extracted from OGR data.
}
  \item{factor.var}{
a boolean vector that indicates variables which have been considered as factor. The length has to be equal to the length of varname argument.
}
}

\details{
The function allows to extract data from OGR data. The argument varname has to contain the name of the ecological response and the names of local variables to be taken into account for the analysis.
}

\value{
a dataframe with location of the observed sites, and the variables (response variable and explanatory variables) used for the analysis.
}

\author{
O. Martin
}

\seealso{
\code{\link[rgdal]{readOGR}}
}


