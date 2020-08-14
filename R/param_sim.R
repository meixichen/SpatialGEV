#' Draw samples of both fixed and random effects in a model with scale parameter as random effect.
#'
#' @param gev_adf A TMB model object that has been optimized 
#' @param y Vector of `n` response values. 
#' @param x `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param n_sim Number of simulations. Default is 1000.
#' @param n_random Number of random effects
#' @param sp.thres Thresholding value to create sparse covariance matrix. Any distance value greater than or equal to `sp.thres` will be set to 0. Default is 0, which means not using sparse matrix.
#' @return A list of two matrices: a matrix of draws from the approximate fixed effect posterior,  
#' and a matrix of draws from the approximate random effect posterior. Each column correspond to one simulation.
#' @export
param_sim <- function(gev_adf, y, x, n_sim = 1000, n_random, sp.thres = 0){
  rep <- TMB::sdreport(gev_adf)
  # get the covariance of fixed effect
  fixed_cov <- rep$cov.fixed 
  # get the estimate of fixed effect
  fixed_rep <- summary(rep, "fixed") 
  fixed_estimate <- fixed_rep[,1]
  
  # Create the data list to be passed on to MakeADFun
  n <- length(y)
  dd <- as.matrix(dist(x))
  data <- list(model="model_b", n=n, y=y, dd=dd, sp_thres=sp.thres)
  
  # preallocate random and fixed effects draws
  random_post_draws <- matrix(rep(NA, n_sim*n_random), nrow = n_random, ncol = n_sim)
  fixed_post_draws <- matrix(rep(NA, n_sim*length(fixed_estimate)), nrow = length(fixed_estimate), ncol = n_sim)
  
  for (i in 1:n_sim){
    flag <- TRUE
    while (flag){
      fixed_draw <- mvtnorm::rmvnorm(1, fixed_estimate, fixed_cov)
      out <- gev_adf$fn(fixed_draw)
      if (!is.na(out)){
        flag <- FALSE
        random_opt <- gev_adf$env$last.par[2:(n_random+1)]
        gev_adf2 <- TMB::MakeADFun(data=data,
                                  parameters=list(a = fixed_draw[1], log_b = random_opt,
                                                  log_s = fixed_draw[2],
                                                  log_sigma = fixed_draw[3], log_ell = fixed_draw[4]),
                                  map = list(a = factor(NA), log_s = factor(NA), log_sigma = factor(NA), log_ell = factor(NA)),
                                  DLL = "SpatialGEV_TMBExports", 
                                  silent = TRUE)
        random_he <- gev_adf2$he(random_opt)
        C <- chol(random_he)
        random_cov <- backsolve(r = C, x = backsolve(r = C, x = diag(nrow(random_he)), transpose = TRUE))
        random_post_draws[, i] <- mvtnorm::rmvnorm(1, random_opt, random_cov)
        fixed_post_draws[, i] <- fixed_draw
      }
    }
  }
  list(fixed = fixed_post_draws,
       random = random_post_draws)
}
