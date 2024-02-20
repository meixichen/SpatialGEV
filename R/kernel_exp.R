#' Exponential covariance function
#'
#' @param x Distance measure.
#' @param sigma The scale parameter with the constraint of `sigma > 0`
#' @param ell The range/lengthscale parameter with the constraint of `ell > 0`.
#' @param X1 A `n1 x 2` matrix containing the coordinates of location set 1.
#' If `x` is not provided, `X1` and `X2` should be provided for calculating their distance.
#' @param X2 A `n2 x 2` coordinate matrix.
#' @return A matrix or a scalar of exponential covariance depending on the type of `x` or
#' whether `X1` and `X2` are used instead.
#' @details Let x = dist(x_i, x_j).
#' ```
#' cov(i,j) = sigma^2*exp(-x/ell)
#' ```
#' @examples
#' X1 <- cbind(runif(10, 1, 10), runif(10, 10, 20))
#' X2 <- cbind(runif(5, 1, 10), runif(5, 10, 20))
#'
#' kernel_exp(sigma=2, ell=1, X1=X1, X2=X2)
#'
#' kernel_exp(as.matrix(stats::dist(X1)), sigma=2, ell=1)
#' @export
kernel_exp <- function(x, sigma, ell, X1=NULL, X2=NULL){
  if (any(c(sigma, ell)<=0)) stop("sigma and ell need to be positive.")
  if (missing(x)){
    if (missing(X1) | missing(X2)) stop("x is not provided. Must provide X1 and X2.")
    if (!is.matrix(X1) | !is.matrix(X2)) stop("X1 and X2 must be matrices.")
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    x <- as.matrix(stats::dist(rbind(X1, X2), method="euclidean"))
    x <- x[1:n1, (n1+1):(n1+n2)]
  }
  sigma^2*exp(-x / ell)
}
