#' Matern covariance function
#'
#' @param x Distance measure.
#' @param phi The range parameter with the constraint of `phi > 0`.
#' @param kappa The smoothness parameter with constraint of `kappa > 0`.
#' @return Matern covariance value.
#' @export
kernel_matern <- function(x, phi, kappa){ 
  if (any(c(phi, kappa)<=0)) stop("phi and kappa need to be positive")
  ifelse(x>0, (2^(kappa-1) * gamma(kappa))^{-1} * (x/phi)^kappa * besselK(x/phi, kappa), 1)
}
