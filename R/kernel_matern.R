#' Matern covariance function
#'
#' @param x Distance measure.
#' @param sigma Positive parameter. (This is in fact sigma^2)
#' @param kappa Positive parameter.
#' @param nu Range parameter default to 1.
#' @return Matern covariance value.
#' @export
kernel_matern <- function(x, sigma, kappa, nu=1){ 
  if (any(c(sigma, kappa)<=0)) stop("sigma and kappa need to be positive")
  ifelse(x>0, 
	 sigma * 2^(1-nu) * (gamma(nu))^{-1} * (kappa*x)^nu * besselK(kappa*x, nu), 
	 sigma)
}
