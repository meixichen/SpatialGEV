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
#' @example examples/r_nll.R
#' @export
r_nll <- function(y, dd, a, log_b, s,
                  hyperparam_a, hyperparam_b, hyperparam_s, kernel = "exp",
                  beta_a = NULL, beta_b = NULL, beta_s = NULL,
                  X_a = NULL, X_b = NULL, X_s = NULL,
                  f_s = function(x) x, ...) {
  n <- length(y)
  if(kernel == "exp") {
    cov_a <- kernel_exp(dd, hyperparam_a[1], hyperparam_a[2], ...)
  } else if(kernel == "matern") {
    cov_a <- kernel_matern(dd, hyperparam_a[1], hyperparam_a[2], ...)
  } else {
    stop("Argument kernel must be `exp` or `matern`.")
  }
  if(is.null(X_a)) X_a <- matrix(1, nrow=n, ncol=1)
  nll <- -mvtnorm::dmvnorm(a, mean = X_a%*%beta_a, sigma = cov_a, log = TRUE)
  # if only a is random
  if(length(log_b)==1 && length(s)==1) {
    for(i in 1:n) {
      nll <- nll - sum(vapply(y[[i]], dgev, loc=a[i], FUN.VALUE=numeric(1),
                              scale=exp(log_b), shape=s, log=TRUE))
    }
  } else if (length(log_b)>1 && length(log_b)==length(a)) {
    # if a and b are both random
    if(missing(hyperparam_b)) {
      stop("b is a spatial random effect. Must provide `hyperparam_b`.")
    }
    if (kernel == "exp") {
      cov_b <- kernel_exp(dd, hyperparam_b[1], hyperparam_b[2], ...)
    } else {
      cov_b <- kernel_matern(dd, hyperparam_b[1], hyperparam_b[2], ...)
    }
    if(is.null(X_b)) X_b <- matrix(1, nrow=n, ncol=1)
    nll <- nll - dmvnorm(log_b, mean = X_b%*%beta_b, sigma = cov_b, log = TRUE)
    if(length(s)==1) {
      # if s is fixed
      for (i in 1:n){
        nll <- nll - sum(vapply(y[[i]], dgev, FUN.VALUE=numeric(1),
				loc=a[i], scale=exp(log_b[i]), shape=s, log=TRUE))
      }
    } else if (length(s)>1 && length(s)==length(a)) {
      if (missing(hyperparam_s)) {
        stop("s is a spatial random effect. Must provide `hyperparam_s`.")
      }
      if(kernel == "exp") {
        cov_s <- kernel_exp(dd, hyperparam_s[1], hyperparam_s[2], ...)
      } else {
        cov_s <- kernel_matern(dd, hyperparam_s[1], hyperparam_s[2], ...)
      }
      if(is.null(X_s)) X_s <- matrix(1, nrow=n, ncol=1)
      nll <- nll - dmvnorm(f_s(s),
                           mean = X_s%*%beta_s, sigma = cov_s, log = TRUE)
      for (i in 1:n){
        nll <- nll - sum(vapply(y[[i]], dgev, FUN.VALUE=numeric(1),
                                loc=a[i], scale=exp(log_b[i]),
                                shape=s[i], log=TRUE))
      }
    } else {
      stop("Check length of `s`: when it is random, it should be same length as `a` and `b`.")
    }
  } else {
    stop("Check the length of `log_b`: when it is random, it should be same length as `a`.")
  }
  nll
}


#' Simulate location data and GP parameters for testing
#' @param random A vector of character strings "a", "b" or "s".
#' @param kernel "exp", "matern", or "spde".
#' @param reparam_s "positive", "negative", "unconstrained", or "zero".
#' @param calc_nll Whether to calculate the negative log-likelihood given the
#' simulated parameters and data? Default to TRUE. If `kernel=="spde"`, this
#' should be FALSE.
#' @return A list of simulated parameters: a, log(b) (if `random` contains "b"),
#' reparameterized s (if `random` contains "s"), corresponding GP hyperparameters,
#' and the negative log-likelihood calculated in R.
test_sim <- function(random="a", kernel=c("exp", "matern", "spde"),
                     reparam_s=c("positive", "negative", "unconstrained", "zero"),
                     calc_nll=TRUE){
  kernel <- match.arg(kernel)
  reparam_s <- match.arg(reparam_s)
  if (kernel == "exp"){
    kernel_fun <- kernel_exp
  } else if (kernel %in% c("matern", "spde")){
    kernel_fun <- kernel_matern
  }
  n_sqrt <- sample(5:10, 1)
  n <- n_sqrt^2
  lon <- seq(0, 10, length.out = n_sqrt)
  lat <- seq(0, 10, length.out = n_sqrt)
  X <- expand.grid(x = lon, y = lat)
  dd <- as.matrix(stats::dist(X))
  gp_hyper1_a <- runif(1, 0, 1)
  gp_hyper2_a <- rnorm(1, 0.5, 0.1)
  cov_a <- kernel_fun(dd, exp(gp_hyper1_a), exp(gp_hyper2_a))
  mean_a <- rep(rnorm(1, 1, 1), n)
  a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
  beta_a <- mean(a)
  f_s <- function(x) x
  if ("b" %in% random){
    gp_hyper1_b <- runif(1, 0,0.5)
    gp_hyper2_b <- rnorm(1, 0.5, 0.1)
    cov_b <- kernel_fun(dd, exp(gp_hyper1_b), exp(gp_hyper2_b))
    mean_b <- rep(rnorm(1, 0.5, 0.5), n)
    log_b <- mvtnorm::rmvnorm(1, mean_b, cov_b)
    beta_b <- mean(log_b)
  } else{
    log_b <- runif(1, -3, 0)
    beta_b <- NULL
    gp_hyper1_b <- NULL
    gp_hyper2_b <- NULL
  }
  if ("s" %in% random){
    gp_hyper1_s <- runif(1, -4, -2)
    gp_hyper2_s <- rnorm(1, 1, 2)
    cov_s <- kernel_fun(dd, exp(gp_hyper1_s), exp(gp_hyper2_s))
    mean_s <- rep(rnorm(1, -3, 1), n)
    log_s <- mvtnorm::rmvnorm(1, mean_s, cov_s)
    s_orig <- exp(log_s)
    beta_s <- mean(log_s)
  } else{
    s_orig <- runif(1, 0.01, 0.1)
    log_s <- log(s_orig)
    tmb_s <- s_orig
    beta_s <- NULL
    gp_hyper1_s <- NULL
    gp_hyper2_s <- NULL
  }
  if (reparam_s == "positive"){
    tmb_s <- log_s
    f_s <- function(x) log(x)
  } else if (reparam_s == "unconstrained"){
    if ("s"%in%random) beta_s <- mean(s_orig)
    tmb_s <- s_orig
  } else if (reparam_s == "negative"){
    s_orig <- -exp(log_s)
    tmb_s <- log_s
    f_s <- function(x) log(abs(s_orig))
  } else if (reparam_s == "zero"){
    s_orig <- 0
    tmb_s <- 0
  }
  y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE),
           loc=a, scale=exp(log_b), shape=s_orig)
  re_list <- list(
    a=a, log_b=log_b, s=tmb_s,
    beta_a=beta_a, beta_b=beta_b, beta_s=beta_s
  )
  hyper_list <- list(
    log_sigma_a=gp_hyper1_a, log_ell_a=gp_hyper2_a,
    log_sigma_b=gp_hyper1_b, log_ell_b=gp_hyper2_b,
    log_sigma_s=gp_hyper1_s, log_ell_s=gp_hyper2_s
  )
  if (kernel %in% c("matern", "spde")){
    names(hyper_list) <- gsub("_ell_", "_kappa_", names(hyper_list))
  }
  param_list <- c(re_list, hyper_list)
  param_list <- param_list[!vapply(param_list, is.null, TRUE)]
  if (calc_nll & kernel!="spde"){
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s_orig,
                   hyperparam_a=c(exp(gp_hyper1_a), exp(gp_hyper2_a)),
                   hyperparam_b=c(exp(gp_hyper1_b), exp(gp_hyper2_b)),
                   hyperparam_s=c(exp(gp_hyper1_s), exp(gp_hyper2_s)),
                   kernel=kernel, beta_a=beta_a, beta_b=beta_b, beta_s=beta_s,
                   f_s=f_s)
  } else{
    nll_r <- NULL
  }
  list(y=y, locs=X, params=param_list, nll_r=nll_r,
       kernel=kernel, random=paste(random, collapse=""), reparam_s=reparam_s)
}


#' Calculate negative log-likelihood in TMB
#' @param sim_res Simulation results produced by `test_sim()`.
#' @return A scalar.
calc_tmb_nll <- function(sim_res){
  adfun <- spatialGEV_fit(sim_res$y, locs=sim_res$locs,
                          random=sim_res$random,
                          init_param=sim_res$params,
                          reparam_s=sim_res$reparam_s,
                          kernel=sim_res$kernel,
                          sp_thres=-1,
                          adfun_only=TRUE,
                          ignore_random=TRUE,
                          silent=TRUE)
  if (sim_res$reparam == "zero"){
    sim_res$params <- sim_res$params[names(sim_res$params) != "s"]
  }
  adfun$fn(unlist(sim_res$params))
}

