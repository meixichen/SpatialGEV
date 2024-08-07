% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simulatedData3}
\alias{simulatedData3}
\title{Simulated dataset 3}
\format{
A list containing the simulation parameters and simulated data on a 20x20 grid:
\describe{
\item{locs}{A 100x2 matrix. First column contains longitudes and second contains latitudes}
\item{a}{A length 100 vector. GEV location parameters}
\item{logb}{A length 100 vector. Log-transformed GEV scale parameters}
\item{logs}{A length 100 vector. Log-transformed GEV shape parameters}
\item{y}{A length 100 list of vectors which are observations simulated at each location}
\item{lon}{A length 10 vector of longitude}
\item{lat}{A length 10 vector of latitude}
}
}
\usage{
simulatedData3
}
\description{
A list of small dataset used for package testing and demos.
\code{a}, \code{logb}, \code{logs} are simulated from smooth surfaces.
}
\keyword{datasets}
