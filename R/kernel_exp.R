#' Exponential covariance function
#'
#' @param x Distance measure.
#' @param sigma The amplitude parameter (scalar) with the constraint of `sigma > 0`
#' @param ell The smoothness parameter (scalar) with the constraint of `ell > 0`.
#' @return Exponential covariance value. 
#' @export
kernel_exp <- function(x, sigma, ell){ 
  if (any(c(sigma, ell)<=0)) stop("sigma and ell need to be positive")
  sigma*exp(-x/ell)
}
