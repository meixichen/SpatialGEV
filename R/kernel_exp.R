#' Exponential covariance function
#'
#' @param X1 A `n1 x d` matrix
#' @param X2 A `n2 x d` matrix
#' @param sigma The amplitude parameter (scalar) with the constraint of `sigma > 0`
#' @param ell The smoothness parameter (scalar) with the constraint of `ell > 0`.
#' @param dist_method A distance measure to be used. See `help(dist)`
#' @return An `n1 x n2` covariance matrix calculated using exponential kernel function.
#' @export
kernel_exp <- function(X1, X2, sigma, ell, dist_method = "euclidean"){ 
  if (any(c(sigma, ell)<=0)) stop("sigma and ell need to be positive")
  if (!is.matrix(X1) | !is.matrix(X2)) stop("x1 and x2 must be matrices")
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dd <- as.matrix(stats::dist(rbind(X1, X2), method = dist_method))
  dd <- dd[1:n1, (n1+1):(n1+n2)]
  sigma*exp(-dd/ell)
}
