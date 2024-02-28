#' Matern covariance function
#'
#' @param x Distance measure.
#' @param sigma Positive scale parameter.
#' @param kappa Positive inverse range/lengthscale parameter.
#' @param nu Smoothness parameter default to 1.
#' @param X1 A `n1 x 2` matrix containing the coordinates of location set 1.
#' If `x` is not provided, `X1` and `X2` should be provided for calculating their distance.
#' @param X2 A `n2 x 2` coordinate matrix.
#' @return A matrix or a scalar of Matern covariance depending on the type of `x` or
#' whether `X1` and `X2` are used instead.
#' @details Let x = dist(x_i, x_j).
#' ```
#' cov(i,j) = sigma^2 * 2^(1-nu)/gamma(nu) * (kappa*x)^nu * K_v(kappa*x)
#' ```
#' Note that when `nu=0.5`, the Matern kernel corresponds to the absolute exponential kernel.
#' @example examples/kernel_matern.R
#' @export
kernel_matern <- function(x, sigma, kappa, nu=1, X1=NULL, X2=NULL){
  if (any(c(sigma, kappa)<=0)) stop("sigma and kappa need to be positive")
  if (missing(x)){
    if (missing(X1) | missing(X2)) stop("x is not provided. Must provide X1 and X2.")
    if (!is.matrix(X1) | !is.matrix(X2)) stop("X1 and X2 must be matrices.")
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    x <- as.matrix(stats::dist(rbind(X1, X2), method="euclidean"))
    x <- x[1:n1, (n1+1):(n1+n2)]
  }
  ifelse(x>0,
	 sigma^2 * 2^(1-nu) * (gamma(nu))^{-1} * (kappa*x)^nu * besselK(kappa*x, nu),
	 sigma^2)
}
