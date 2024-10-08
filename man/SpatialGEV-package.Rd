% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialGEV-package.R
\docType{package}
\name{SpatialGEV-package}
\alias{SpatialGEV}
\alias{SpatialGEV-package}
\title{SpatialGEV: Fit Spatial Generalized Extreme Value Models}
\description{
Fit latent variable models with the GEV distribution as the data likelihood and the GEV parameters following latent Gaussian processes. The models in this package are built using the template model builder 'TMB' in R, which has the fast ability to integrate out the latent variables using Laplace approximation. This package allows the users to choose in the fit function which GEV parameter(s) is considered as a spatially varying random effect following a Gaussian process, so the users can fit spatial GEV models with different complexities to their dataset without having to write the models in 'TMB' by themselves. This package also offers methods to sample from both fixed and random effects posteriors as well as the posterior predictive distributions at different spatial locations. Methods for fitting this class of models are described in Chen, Ramezan, and Lysy (2024) \doi{10.1002/cjs.11829}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/meixichen/SpatialGEV}
  \item Report bugs at \url{https://github.com/meixichen/SpatialGEV/issues}
}

}
\author{
\strong{Maintainer}: Meixi Chen \email{meixi.chen@uwaterloo.ca}

Authors:
\itemize{
  \item Martin Lysy \email{mlysy@uwaterloo.ca}
}

Other contributors:
\itemize{
  \item Reza Ramezan \email{rramezan@uwaterloo.ca} [contributor]
}

}
\keyword{internal}
