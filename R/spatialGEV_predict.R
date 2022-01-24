#' Draw from the posterior predictive distributions at new locations based on a fitted spatialGEV model
#'
#' @param model A fitted spatial GEV model returned by `spatialGEV_fit`
#' @param X_new A `n_test x 2` matrix containing the coordinates of the new locations
#' @param X_obs A `n_obs x 2` matrix containing the coordinates of the observed locations for model fitting
#' @param n_draw Number of draws from the posterior predictive distribution
#' @return An `n_draw x n_test` matrix containing the draws from the posterior predictive distributions at `n_test` new locations
#' @export
spatialGEV_predict <- function(model, X_new, X_obs, n_draw){
  # extract info from model
  n_test <- nrow(X_new)
  rep <- model$report
  joint_mean <- c(rep$par.random, rep$par.fixed)
  random_ind <- rep$env$random #indices of random effects
  fixed_ind <- (1:length(joint_mean))[-random_ind] #indices of fixed effects
  n_obs <- nrow(model$adfun$env$data$dd) # number of locations
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  if (length(random_ind) == n_obs) {
    mod <- "a"
  }
  else if (length(random_ind) == 2*n_obs){
    mod <- "ab"
  }
  else {
    stop("n_obs must divide the length of random effect vector.")
  }
  # Sample from MVN
  jointPrec_mat <- rep$jointPrecision
  C <- chol(jointPrec_mat)
  joint_cov <- backsolve(r = C, x = backsolve(r = C, x = diag(nrow(jointPrec_mat)), transpose = TRUE))
  parameter_draw <- mvtnorm::rmvnorm(n_draw, joint_mean, joint_cov)
  
  s_draw_fun <- function(j, s_ind){ 
    # `j` points to the j-th draw and `s_ind` is the index of s in the parameter vector of the model
    if (reparam_s == 3){
      parameter_draw[j, s_ind]
    }
    else if (reparam_s == 1){
      exp(parameter_draw[j, s_ind])
    }
    else if (reparam_s == 2){
      -exp(parameter_draw[j, s_ind])
    }
    else{
      0
    }
  }
  total_param <- ncol(parameter_draw)
  pred_y_draws <- matrix(NA, nrow = n_draw, ncol = n_test)
  if (mod == "a"){
    for (i in 1:n_draw){
      a <- parameter_draw[i, 1:n_obs]
      b <- exp(parameter_draw[i, n_obs+1])
      s <- s_draw_fun(i, n_obs+2)
      sigma_a <- exp(parameter_draw[i, total_param-1])
      ell_a <- exp(parameter_draw[i, total_param])
      # Construct conditional distribution function for a
      a_sim_fun <- sim_cond_normal(rep(0, (n_obs+n_test)), a = a, X.new = X_new, X.obs = as.matrix(X_obs), 
                                   kernel = kernel_exp, sigma = sigma_a, ell = ell_a)
      new_a <- a_sim_fun(1)
      new_y <- t(apply(X = new_a, MARGIN = 1, FUN = function(row){
        unlist(Map(rgev, n=1, loc=row, scale=b, shape=s))  
      })) # a `1 x n_test` matrix
      pred_y_draws[i, ] <- new_y
    }
  }
  else {
    for (i in 1:n_draw){
      a <- parameter_draw[i, 1:n_obs]
      logb <- parameter_draw[i, (n_obs+1):(2*n_obs)]
      s <- s_draw_fun(i, 2*n_obs + 1)
      sigma_a <- exp(parameter_draw[i, total_param-3])
      ell_a <- exp(parameter_draw[i, total_param-2])
      sigma_b <- exp(parameter_draw[i, total_param-1])
      ell_b <- exp(parameter_draw[i, total_param])
      # Construct conditional distribution function for a and logb
      a_sim_fun <- sim_cond_normal(rep(0, (n_obs+n_test)), a = a, X.new = as.matrix(X_new), X.obs = as.matrix(X_obs), 
                                   kernel = kernel_exp, sigma = sigma_a, ell = ell_a)
      logb_sim_fun <- sim_cond_normal(rep(0, (n_obs+n_test)), a = logb, X.new = as.matrix(X_new), X.obs = as.matrix(X_obs), 
                                      kernel = kernel_exp, sigma = sigma_b, ell = ell_b)
      
      new_a <- a_sim_fun(1) # 1 x n_test matrix
      new_logb <- logb_sim_fun(1) # 1 x n_test matrix
      new_ab <- cbind(new_a, exp(new_logb)) # A `1 x (2*n_test)` matrix constructed by putting the matrix of exp(logb) to the right of the matrix of a
      new_y <- t(apply(X = new_ab, MARGIN = 1, FUN = function(row){
        unlist(Map(rgev, n=1, loc=row[1:n_test], scale=row[(n_test+1):length(row)], shape=s))  
      })) # a `1 x n_test` matrix
      pred_y_draws[i, ] <- new_y
    }
  }
  pred_y_draws
}
