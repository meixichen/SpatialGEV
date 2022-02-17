#' Matern covariance function
#'
#' @param x Distance measure.
#' @param sigma Positive parameter. (This is in fact sigma^2)
#' @param kappa Positive parameter.
#' @param nu Range parameter default to 1.
#' @param X1 A `n1 x 2` matrix containing the coordinates of location set 1. 
#' If `x` is not provided, `X1` and `X2` should be provided for calculating their distance.
#' @param X2 A `n2 x 2` coordinate matrix.
#' @return A matrix or a scalar of Matern covariance depending on the type of `x` or 
#' whether `X1` and `X2` are used instead. 
#' @details Let x = dist(x_i, x_j).
#' ```
#' cov(i,j) = sigma * 2^(1-nu)/gamma(nu) * (kappa*x)^nu * K_v(kappa*x)
#' ```
#' @examples
#' X1 <- cbind(runif(10, 1, 10), runif(10, 10, 20))
#' X2 <- cbind(runif(5, 1, 10), runif(5, 10, 20))
#'
#' kernel_matern(sigma=2, ell=1, X1=X1, X2=X2)
#'
#' kernel_matern(as.matrix(stats::dist(X1)), sigma=2, ell=1)
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
	 sigma * 2^(1-nu) * (gamma(nu))^{-1} * (kappa*x)^nu * besselK(kappa*x, nu), 
	 sigma)
}
