#' Calculate the negative marginal loglikelihood of the GEV-GP  model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations.
#' @param dd An `n x n` distance matrix.
#' @param a Vector of `n` location paramter
#' @param log_b A numeric value or a vector of `n` log-transformed scale parameters if considered 
#' as a random effect.
#' @param s A numeric value or a vector of `n` shape parameters
#' @param hyperparam_a A vector of hyperparameters for a. See details.
#' @param hyperparam_b A vector of hyperparameters for b. Must be provided if `log_b` is a vector. 
#' See details.
#' @param hyperparam_s A vector of hyperparameters for f(s), where f() is a transformation function
#' for s specifided using the `f_s` argument. Must be provided if `s` is a vector. 
#' @param kernel "exp" or "matern". Kernel function used to compute the covariance matrix for 
#' spatial random effects. Default is "exp".
#' @param beta_a Numeric. Coefficients for mean of GP(a).
#' @param beta_b Numeric. Coefficients for mean of GP(log_b).
#' @param beta_s Numeric. Coefficients for mean of GP(s).
#' @param X_a Design matrix for a. If not provided, this will a `n_loc x 1` column matrix of 1s.
#' @param X_b Design matrix for log(b). If not provided and logb is a random effect, 
#' this will a `n_loc x 1` column matrix of 1s.
#' @param X_s Design matrix for s. If not provided, this will a `n_loc x 1` column matrix of 1s.
#' @param f_s A function f() used to transform s such that 
#' f(s) ~ GP(X_s*beta_s, Sigma(hyperparam_s)). Default is identitfy function: `function(x){x}`.
#' @param ... Additional arguments to pass to the kernel function, e.g. `nu` for the matern.
#' @return Scalar value of the negative marginal loglikelihood:
#' ```
#' -logL(Data; spatial_random_effects, fixed_hyperparameters)
#' ```
#' @details This function is used to test if TMB and R output the same negative loglikelihood.
#' If `kernel="exp`, `hyperparam_a/b/s` should be `c(sigma_a/b/s, ell_a/b/s)`, where `sigma` is the 
#' amplitude hyperparameter and `ell` is the smoothness hyperparameter for the exponential kernel.
#' If `kernel="matern`, `hyperparam_a/b/s` should be `c(sigma_a/b, kappa_a/b/s)`, where `sigma` and 
#' `kappa` are hyperparameters for the Matern kernel. 
#' If only `a` is a spatial random effect and `b` is fixed, only `hyperparam_a` needs to be 
#' provided.
#'
#' This function is used as the ground truth for testing hpp model likelihood.
#' @examples
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' s <- exp(logs)
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' dd <- as.matrix(stats::dist(locs))
#' log_sigma_a <- -1; log_ell_a <- 5
#' log_sigma_b <- -2; log_ell_b <- 10
#' beta_a <- mean(a); beta_b <- mean(logb)
#' # Negative marginal log-likelihood produced in R using the exponential kernel
#' nll_r <- r_nll(y, dd, a=a, log_b=logb, s=s,
#'                hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
#'                hyperparam_b=c(exp(log_sigma_b), exp(log_ell_b)),
#'                kernel="exp", beta_a=beta_a, beta_b=beta_b)
#' # Negative marg loglik produced by TMB template
#' init_param <- list(beta_a=beta_a, beta_b=beta_b,
#'                    a=a, log_b=logb, s=log(s), 
#'                    log_sigma_a=log_sigma_a, 
#'                    log_ell_a=log_ell_a,
#'                    log_sigma_b=log_sigma_b, 
#'                    log_ell_b=log_ell_b)
#' adfun <- spatialGEV_fit(y, locs, random="ab",
#'                         init_param=init_param,
#'                         reparam_s="positive",
#'                         kernel="exp",
#'                         adfun_only=TRUE,
#'                         ignore_random=TRUE,
#'                         silent=TRUE)
#' nll_tmb <- adfun$fn(unlist(init_param))
#' nll_r - nll_tmb
#' @export
r_nll <- function(y, dd, a, log_b, s, hyperparam_a, hyperparam_b, hyperparam_s, kernel = "exp", 
		  beta_a = NULL, beta_b = NULL, beta_s = NULL, X_a = NULL, X_b = NULL, X_s = NULL, 
		  f_s = function(x){x}, ...) {
  n <- length(y)
  if (kernel == "exp"){
    cov_a <- kernel_exp(dd, hyperparam_a[1], hyperparam_a[2], ...)
  }
  else if (kernel == "matern"){
    cov_a <- kernel_matern(dd, hyperparam_a[1], hyperparam_a[2], ...)
  }
  else{
    stop("Argument kernel must be `exp` or `matern`.")
  }
  if (is.null(X_a)) X_a <- matrix(1, nrow=n, ncol=1)
  nll <- -mvtnorm::dmvnorm(a, mean = X_a%*%beta_a, sigma = cov_a, log = TRUE) 
  # if only a is random
  if (length(log_b)==1 & length(s)==1){ 
    for (i in 1:n){
      nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b), shape=s, log=TRUE))
    }
  }
  # if a and b are both random
  else if (length(log_b)>1 & length(log_b)==length(a)){
    if (missing(hyperparam_b)){stop("b is a spatial random effect. Must provide `hyperparam_b`.")}
    if (kernel == "exp"){
      cov_b <- kernel_exp(dd, hyperparam_b[1], hyperparam_b[2], ...)
    }else{
      cov_b <- kernel_matern(dd, hyperparam_b[1], hyperparam_b[2], ...)
    }
    if (is.null(X_b)) X_b <- matrix(1, nrow=n, ncol=1)
    nll <- nll - dmvnorm(log_b, mean = X_b%*%beta_b, sigma = cov_b, log = TRUE)
    
    # if s is fixed
    if (length(s)==1){
      for (i in 1:n){
	nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b[i]), shape=s, log=TRUE))
      }
    }
    else if (length(s)>1 & length(s)==length(a)){
      if (missing(hyperparam_s)){stop("s is a spatial random effect. Must provide `hyperparam_s`.")}
      if (kernel == "exp"){
	cov_s <- kernel_exp(dd, hyperparam_s[1], hyperparam_s[2], ...)
      }else{
	cov_s <- kernel_matern(dd, hyperparam_s[1], hyperparam_s[2], ...)
      }
      if (is.null(X_s)) X_s <- matrix(1, nrow=n, ncol=1)
      nll <- nll - dmvnorm(f_s(s), mean = X_s%*%beta_s, sigma = cov_s, log = TRUE)
      for (i in 1:n){
	nll <- nll - sum(sapply(y[[i]], dgev, 
				loc=a[i], scale=exp(log_b[i]), shape=s[i], log=TRUE))
      }
    }
    else{
      stop("Check length of `s`: when it is random, it should be same length as `a` and `b`.")
    }
  }
  else{
    stop("Check the length of `log_b`: when it is random, it should be same length as `a`.")
  }

  nll
}
