#' Calculate the negative loglikelihood of the spatial GEV model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations.
#' @param dd An `n x n` distance matrix.
#' @param a Vector of `n` location paramter
#' @param log_b A numeric value or a vector of `n` log-transformed scale parameters if considered as a random effect.
#' @param s Shape parameter
#' @param hyperparam_a A vector of hyperparameters for a. See details.
#' @param hyperparam_b A vector of hyperparameters for b. Must be provided if `log_b` is a vector. See details.
#' @param kernel "exp" or "matern". Kernel function used to compute the covariance matrix for spatial random effects. Default is "exp". 
#' @return Scalar value of the negative loglikelihood.
#' @details This function is used to test if TMB and R output the same negative loglikelihood.
#' If `kernel="exp`, `hyperparam_a/b` should be `c(sigma_a/b, ell_a/b)`, where `sigma` is the amplitude hyperparameter and `ell` is the smoothness hyperparameter for the exponential kernel.
#' If `kernel="matern`, `hyperparam_a/b` should be `c(phi_a/b, kappa_a/b)`, where `phi` is the range hyperparameter and `kappa` is the smoothness hyperparameter for the Matern kernel. 
#' If only `a` is a spatial random effect and `b` is fixed, only `hyperparam_a` needs to be provided.
#' @export
r_nll <- function(y, dd, a, log_b, s, hyperparam_a, hyperparam_b, kernel = "exp") {
  n <- length(y)
  if (kernel == "exp"){
    cov_a <- kernel_exp(dd, hyperparam_a[1], hyperparam_a[2])
  }else if (kernel == "matern"){
    cov_a <- kernel_matern(dd, hyperparam_a[1], hyperparam_a[2])
  }else{
    stop("Argument kernel must be `exp` or `matern`.")
  }
  nll <- -dmvnorm(a, sigma = cov_a, log = TRUE) 
  if (length(log_b)==1){
    for (i in 1:n){
      nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b), shape=s, log=TRUE))
    }
  }
  else { # if b is considered a random effect
    if (missing(hyperparam_b)){stop("b is a spatial random effect. Must provide `hyperparam_b`.")}
    if (kernel == "exp"){
      cov_b <- kernel_exp(dd, hyperparam_b[1], hyperparam_b[2])
    }else{
      cov_b <- kernel_matern(dd, hyperparam_b[1], hyperparam_b[2])
    }
    nll <- nll - dmvnorm(log_b, sigma = cov_b, log = TRUE)
    for (i in 1:n){
      nll <- nll - sum(sapply(y[[i]], dgev, loc=a[i], scale=exp(log_b[i]), shape=s, log=TRUE))
    }
  }
  nll
}
