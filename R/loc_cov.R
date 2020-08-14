#' The spatial exponential covariance function
#'
#' @param dd An `n x n` matrix for the distances between `n` locations.
#' @param sigma The amplitude parameter (scalar) with the constraint of `sigma > 0`
#' @param ell The smoothness parameter (scalar) with the constraint of `ell > 0`.
#' @return An `n x n` covariance matrix.
#' @export
loc_cov <- function(dd, sigma, ell){ 
  if (any(c(sigma, ell)<=0)) stop("sigma and ell need to be positive")
  n <- nrow(dd)
  sigma*exp(-dd/ell)
}
