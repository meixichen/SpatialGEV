#' Calculate the negative marginal loglikelihood of the GEV-GP  model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations.
#' @param dd An `n x n` distance matrix.
#' @param a Vector of `n` location paramter
#' @param log_b A numeric value or a vector of `n` log-transformed scale parameters if considered 
#' as a random effect.
#' @param s Shape parameter
#' @param hyperparam_a A vector of hyperparameters for a. See details.
#' @param hyperparam_b A vector of hyperparameters for b. Must be provided if `log_b` is a vector. 
#' See details.
#' @param kernel "exp" or "matern". Kernel function used to compute the covariance matrix for 
#' spatial random effects. Default is "exp".
#' @param beta_a Numeric. Coefficients for mean of GP(a).
#' @param beta_b Numeric. Coefficients for mean of GP(log_b).
#' @param X_a Design matrix for a. If not provided, this will a `n_loc x 1` column matrix of 1s.
#' @param X_b Design matrix for log(b). If not provided and logb is a random effect, 
#' this will a `n_loc x 1` column matrix of 1s.
#' @param ... Additional arguments to pass to the kernel function, e.g. `nu` for the matern.
#' @return Scalar value of the negative marginal loglikelihood:
#' ```
#' -logL(Data; spatial_random_effects, fixed_hyperparameters)
#' ```
#' @details This function is used to test if TMB and R output the same negative loglikelihood.
#' If `kernel="exp`, `hyperparam_a/b` should be `c(sigma_a/b, ell_a/b)`, where `sigma` is the 
#' amplitude hyperparameter and `ell` is the smoothness hyperparameter for the exponential kernel.
#' If `kernel="matern`, `hyperparam_a/b` should be `c(sigma_a/b, kappa_a/b)`, where `sigma` and 
#' `kappa` are hyperparameters for the Matern kernel. 
#' If only `a` is a spatial random effect and `b` is fixed, only `hyperparam_a` needs to be 
#' provided.
#' @examples
#' \dontrun{
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' n_loc <- nrow(locs)
#' log_sigma_a <- -1; log_ell_a <- 5
#' log_sigma_b <- -2; log_ell_b <- 10
#' beta_a <- mean(a); beta_b <- mean(logb)
#' # Negative marginal log-likelihood produced in R using the exponential kernel
#' nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
#'                hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
#'                hyperparam_b=c(exo(log_sigma_b), exp(log_ell_b)),
#'                kernel="exp", beta_a=beta_a, beta_b=beta_b)
#' # Negative marg loglik produced by TMB template
#' adfun <- spatialGEV_fit(y, X, random="ab",
#'                         init_param=list(beta_a=beta_a, beta_b=beta_b,
#'                                         a=a, log_b=log_b, s=log(s), 
#'                                         log_sigma_a=log_sigma_a, 
#'                                         log_ell_a=log_ell_a,
#'                                         log_sigma_b=log_sigma_b, 
#'                                         log_ell_b=log_ell_b)
#'                         reparam_s="positive",
#'                         kernel="exp",
#'                         adfun_only=TRUE,
#'                         ignore_random=TRUE,
#'                         silent=TRUE)
#' nll_tmb <- adfun$fn(unlist(init_param))
#' testthat::expect_equal(nll_r, nll_tmb)
#' }
#' @export
r_nll <- function(y, dd, a, log_b, s, hyperparam_a, hyperparam_b, kernel = "exp", 
		  beta_a = NULL, beta_b = NULL, X_a = NULL, X_b = NULL, ...) {
  n <- length(y)
  if (kernel == "exp"){
    cov_a <- kernel_exp(dd, hyperparam_a[1], hyperparam_a[2], ...)
  }else if (kernel == "matern"){
    cov_a <- kernel_matern(dd, hyperparam_a[1], hyperparam_a[2], ...)
  }else{
    stop("Argument kernel must be `exp` or `matern`.")
  }
  if (is.null(X_a)) X_a <- matrix(1, nrow=n, ncol=1)
  nll <- -mvtnorm::dmvnorm(a, mean = X_a%*%beta_a, sigma = cov_a, log = TRUE) 
  if (length(log_b)==1){
    for (i in 1:n){
      nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b), shape=s, log=TRUE))
    }
  }
  else { # if b is considered a random effect
    if (missing(hyperparam_b)){stop("b is a spatial random effect. Must provide `hyperparam_b`.")}
    if (kernel == "exp"){
      cov_b <- kernel_exp(dd, hyperparam_b[1], hyperparam_b[2], ...)
    }else{
      cov_b <- kernel_matern(dd, hyperparam_b[1], hyperparam_b[2], ...)
    }
    if (is.null(X_b)) X_b <- matrix(1, nrow=n, ncol=1)
    nll <- nll - dmvnorm(log_b, mean = X_b%*%beta_b, sigma = cov_b, log = TRUE)
    for (i in 1:n){
      nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b[i]), shape=s, log=TRUE))
    }
  }
  nll
}
