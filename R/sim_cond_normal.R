#' Create a function to simulate from the conditional normal distribution of new data given old data
#'
#' @param joint.mean The length `n` mean vector of the MVN distribution. By default mu1 is the first `m` elements of `joint.mean`
#' @param a A vector of length `n-m`, the values of mu2 to condition on
#' @param X.new A matrix containing the coordiantes of new locations
#' @param X.obs A matrix containing the coordinates of observed locations
#' @param kernel A function that returns a matrix containing the similarity between the two arguments 
#' @param sigma Parameter of kernel function
#' @param ell Parameter of kernel function
#' @return A function that takes in one argument `n` as the number of samples to draw from the condition normal distribution
#' of `X.new` given `X.old`: either from `rmvnorm` for MVN or `rnorm` for univariate normal. The old and new data are assumed
#' to follow a joint multivariate normal distribution. 
#' @details The notations are consistent to the notations on the MVN wikipedia page
#' @export
sim_cond_normal <- function(joint.mean, a, X.new, X.obs, kernel=kernel_exp, sigma, ell){
  if (!is.matrix(X.new) | !is.matrix(X.obs)) stop("X.new and X.obs must be matrices")
  n <- length(joint.mean)
  m <- nrow(X.new)
  if (m < 1 | m >= n) stop("Invalid length of mu1")
  mu1 <- joint.mean[1:m]
  mu2 <- joint.mean[(m+1):n]
  Sig11 <- kernel(X1 = X.new, X2 = X.new, sigma = sigma, ell = ell)
  Sig12 <- kernel(X1 = X.new, X2 = X.obs, sigma = sigma, ell = ell)
  if (!is.matrix(Sig12)) Sig12 <- matrix(Sig12, nrow = m)
  Sig21 <- t(Sig12)
  Sig22 <- kernel(X1 = X.obs, X2 = X.obs, sigma = sigma, ell = ell)
  
  # Matrix inversion using cholesky decomposition
  C <- chol(Sig22)
  A.t <- backsolve(r = C, x = backsolve(r = C, x = Sig21, transpose = TRUE)) # A.t = Sig22^{-1} * Sig21
  A <- t(A.t) #A = Sig12 * Sig22^{-1}
  
  mu.bar <- mu1 + A %*% (a - mu2)
  Sig.bar <- Sig11 - Sig12 %*% A.t
  
  if (m == 1) {
    return( function(n) {rnorm(n, mean = mu.bar, sd = sqrt(Sig.bar))} )
  }
  else {
    Sig.bar[lower.tri(Sig.bar)] <- t(Sig.bar)[lower.tri(Sig.bar)] # map the upper triangle to the lower triangle in case Sig.bar is not symmetric due to numerical cancellation
    return( function(n) {mvtnorm::rmvnorm(n, mean = mu.bar, sigma = Sig.bar)})
  }
}
