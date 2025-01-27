% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spilhaus.R
\name{spilhaus_lonlat}
\alias{spilhaus_lonlat}
\title{Project points from Spilhaus to longlat}
\usage{
spilhaus_lonlat(x)
}
\arguments{
\item{x}{matrix of coordinates in Spilhaus}
}
\value{
matrix of longlat coordinates
}
\description{
Input a matrix of Spilhaus coordinates and return longlat coordinates.
}
\details{
Note that this transformation is not available via standard libraries, and was
implemented in R and Python by Ricardo T. Lemos.
}
\examples{
spilhaus_lonlat(cbind(0, 0))
}
