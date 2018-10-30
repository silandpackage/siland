\name{land.gis}
\alias{land.gis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read OGR vector maps for landscape descritpion.
}
\description{
The function reads an OGR data source and return a list of dataframes. The function is developped to 
read easily landscape description from OGR data source. A similar function (data.gis) 
is available to read OGR data source for oberved data or measurements for an ecological variable.


}
\usage{
land.gis(dsn, layer, varname, landname, wd = 100, extentLand = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dsn}{
data source name (interpretation varies by driver, for some drivers, dsn is a file name, but may also be a folder)
}
  \item{layer}{
layer name (varies by driver, may be a file name without extension). From rgdal 1.2.*, layer may be missing, in which case ogrListLayers examines the dsn, and fails if there are no layers, silently reads the only layer if only one layer is found, and reads the first layer if multiple layers are present, issuing a warning that layer should be given explicitly.
}
  \item{varname}{
name of the variable that contains the landscape features
}
  \item{landname}{
a vector that contains the landscape features that are necessary to extract from the OGR vector maps
}
  \item{wd}{
the size of pixels for x-axis and y-axis. By default, wd=100.
}
  \item{extentLand}{
an extent object. If argument extentLand is NULL (by default), the extent from OGR data is retained. In 
other case, the landscape description is limited to the extent given in this argument.

}
}

\value{
The returned value is a list of dataframes. The first element of 
the list corresponds to the the first landscape feature of the argument landname, and so on. In each 
dataframe, pixel locations associated to a landscape feature are given.
}

\author{
O. Martin
}


\section{Warning }{The function can generate some warnings due to the use of NA values to build outputs lists. Do not take into accont warnings concerning these warnings.}

\seealso{
\code{\link[rgdal]{readOGR}}, \code{\link[raster]{extent}}

}

