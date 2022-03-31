#' Draw from the posterior predictive distributions at new locations based on a fitted GEV-GP model
#'
#' @param model A fitted spatial GEV model object of class `spatialGEVfit`
#' @param locs_new A `n_test x 2` matrix containing the coordinates of the new locations
#' @param n_draw Number of draws from the posterior predictive distribution
#' @param X_a_new `n_test x r1` design matrix for a at the new locations. If not provided, the 
#' default is a column matrix of all 1s.
#' @param X_b_new `n_test x r2` design matrix for log(b) at the new locations 
#' @param X_s_new `n_test x r2` design matrix for (possibly transformed) s at the new locations 
#' @param parameter_draws Optional. A `n_draw x n_parameter` matrix. If `spatialGEV_sample()` has 
#' already been called, the output matrix of parameter draws can be supplied here to avoid doing
#' sampling of parameters again. Make sure the number of rows of `parameter_draws` is the same as
#' `n_draw`.
#' @return An object of class `spatialGEVpred`, which is a list of the following components: 
#' - An `n_draw x n_test` matrix `pred_y_draws` containing the draws from the posterior predictive 
#' distributions at `n_test` new locations
#' - An `n_test x 2` matrix `locs_new` containing the coordinates of the test data
#' - An `n_train x 2` matrix `locs_obs` containing the coordinates of the observed data
#' @examples 
#' \donttest{
#' set.seed(123)
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' n_loc <- nrow(locs)
#' n_test <- 20
#' test_ind <- sample(1:n_loc, n_test)
#'
#' # Obtain coordinate matrices and data lists
#' locs_test <- locs[test_ind,]
#' y_test <- y[test_ind]
#' locs_train <- locs[-test_ind,]
#' y_train <- y[-test_ind]
#'
#' # Fit the GEV-GP model to the training set
#' train_fit <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
#'                             init_param = list(beta_a = mean(a),
#'                                              beta_b = mean(logb),   
#'                                              a = rep(0, n_loc-n_test), 
#'	              				log_b = rep(0, n_loc-n_test),
#'						s = 0,
#'						log_sigma_a = 1, 
#'                                              log_kappa_a = -2,
#'						log_sigma_b = 1, 
#'                                              log_kappa_b = -2),
#'		       	       reparam_s = "positive", 
#'			       kernel = "matern",
#'			       silent = TRUE)
#' pred <- spatialGEV_predict(model = train_fit, locs_new = locs_test, n_draw = 100)
#' summary(pred)
#' }
#' @export
spatialGEV_predict <- function(model, locs_new, n_draw, X_a_new=NULL, X_b_new=NULL, X_s_new=NULL, 
			       parameter_draws=NULL){
  # extract info from model
  locs_obs <- model$locs_obs
  X_a <- model$X_a
  X_b <- model$X_b
  X_s <- model$X_s
  kernel <- model$kernel
  nu <- model$nu # Matern hyperparameter
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  n_test <- nrow(locs_new)
  n_train <- length(model$adfun$env$data$n_obs) 
  random_ind <- model$rep$env$random # indices of random effects
  random <- model$random
  if (length(random) == 1) {
    mod <- "a"
  }
  else if (length(random) == 2){
    mod <- "ab"
  }
  else if (length(random) == 3){
    mod <- "abs"
  }
  else {
    stop("Cannot identify which GEV parameters are random.")
  }


  # get parameter draws
  if (is.null(parameter_draws)){ 
    parameter_draws <- spatialGEV_sample(model, n_draw, observation=FALSE)$parameter_draws
  }
  else{
    if (nrow(parameter_draws) != n_draw) stop("nrow(parameter_draws) must be n_draw.")
  }

  s_draw_fun <- function(j, s_ind){ 
    # `j` points to the j-th draw and `s_ind` is the index of s in the parameter vector of the model
    if (reparam_s == 3){
      parameter_draws[j, s_ind]
    }
    else if (reparam_s == 1){
      exp(parameter_draws[j, s_ind])
    }
    else if (reparam_s == 2){
      -exp(parameter_draws[j, s_ind])
    }
    else{
      0
    }
  }
  total_param <- ncol(parameter_draws)
  pred_y_draws <- matrix(NA, nrow = n_draw, ncol = n_test)
  s_ind <- which(colnames(parameter_draws)=="s")

  # Sampling depends on model type
  if (is.null(X_a_new)) X_a_new <- matrix(1, nrow=n_test, ncol=1) # Default design matrix for a
  if (ncol(X_a_new) != ncol(X_a)) stop("Dimensions of X_a_new and X_a must match.") 
  beta_a_ind <- grep("beta_a", colnames(parameter_draws)) # indices of betas for a
  if (mod == "a"){
    for (i in 1:n_draw){
      a <- parameter_draws[i, 1:n_train]
      b <- exp(parameter_draws[i, n_train+1])
      beta_a <- parameter_draws[i, beta_a_ind]
      s <- s_draw_fun(i, s_ind)
      X_all <- rbind(X_a, X_a_new)
      # Construct conditional distribution function for a
      if (kernel == "exp"){
        a_sim_fun <- sim_cond_normal(X_all%*%beta_a, a = a, 
	  			     locs_new = locs_new, locs_obs = as.matrix(locs_obs), 
                                     kernel = kernel_exp, 
				     sigma = exp(parameter_draws[i,"log_sigma_a"]), 
				     ell = exp(parameter_draws[i,"log_ell_a"]))
      }
      else{ 
        if (kernel == "spde"){
	  meshidxloc <- model$meshidxloc
          X_all <- rbind(as.matrix(X_a[meshidxloc,]), X_a_new)
	}
	a_sim_fun <- sim_cond_normal(X_all%*%beta_a, a = a,
                                     locs_new = locs_new, locs_obs = as.matrix(locs_obs),
                                     kernel = kernel_matern, 
				     sigma = exp(parameter_draws[i,"log_sigma_a"]) , 
				     kappa = exp(parameter_draws[i,"log_kappa_a"]), 
				     nu = nu)
      }
      new_a <- a_sim_fun(1) # sample parameter a one time
      new_y <- t(apply(X = new_a, MARGIN = 1, FUN = function(row){
        unlist(Map(rgev, n=1, loc=row, scale=b, shape=s))  
      })) # a `1 x n_test` matrix
      pred_y_draws[i, ] <- new_y
    }
  }
  else if (mod == "ab"){ 
    if (is.null(X_b_new)) X_b_new <- matrix(1, nrow=n_test, ncol=1) # Default design matrix for a
    if (ncol(X_b_new) != ncol(X_b)) stop("Dimensions of X_b_new and X_b must match.") 
    beta_b_ind <- grep("beta_b", colnames(parameter_draws)) # indices of betas for b
    for (i in 1:n_draw){
      a <- parameter_draws[i, 1:n_train]
      logb <- parameter_draws[i, (n_train+1):(2*n_train)]
      beta_a <- parameter_draws[i, beta_a_ind]
      beta_b <- parameter_draws[i, beta_b_ind]
      s <- s_draw_fun(i, s_ind)
      X_all_a <- rbind(X_a, X_a_new)
      X_all_b <- rbind(X_b, X_b_new)
      hyperparam_a1 <- exp(parameter_draws[i, "log_sigma_a"])
      hyperparam_b1 <- exp(parameter_draws[i, "log_sigma_b"])
      # Construct conditional distribution function for a and logb
      if (kernel == "exp"){
	hyperparam_a2 <- exp(parameter_draws[i, "log_ell_a"])
	hyperparam_b2 <- exp(parameter_draws[i, "log_ell_b"])
	a_sim_fun <- sim_cond_normal(X_all_a%*%beta_a, a = a, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_exp, 
				     sigma = hyperparam_a1, 
				     ell = hyperparam_a2)
	logb_sim_fun <- sim_cond_normal(X_all_b%*%beta_b, a = logb, 
					locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
					kernel = kernel_exp, 
					sigma = hyperparam_b1, 
					ell = hyperparam_b2)
      }
      else{
	hyperparam_a2 <- exp(parameter_draws[i, "log_kappa_a"])
	hyperparam_b2 <- exp(parameter_draws[i, "log_kappa_b"])
        if (kernel == "spde"){
	  meshidxloc <- model$meshidxloc
          X_all_a <- rbind(as.matrix(X_a[meshidxloc,]), X_a_new)
          X_all_b <- rbind(as.matrix(X_b[meshidxloc,]), X_b_new)
	}
	a_sim_fun <- sim_cond_normal(X_all_a%*%beta_a, a = a, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_matern, sigma = hyperparam_a1, 
				     kappa = hyperparam_a2, nu = nu)
	logb_sim_fun <- sim_cond_normal(X_all_b%*%beta_b, a = logb, 
					locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
					kernel = kernel_matern, sigma = hyperparam_b1, 
					kappa = hyperparam_b2, nu = nu)
      }
      new_a <- a_sim_fun(1) # 1 x n_test matrix
      new_logb <- logb_sim_fun(1) # 1 x n_test matrix
      new_ab <- cbind(new_a, exp(new_logb)) # A `1 x (2*n_test)` matrix constructed by putting the matrix of exp(logb) to the right of the matrix of a
      new_y <- t(apply(X = new_ab, MARGIN = 1, FUN = function(row){
        unlist(Map(rgev, n=1, loc=row[1:n_test], scale=row[(n_test+1):length(row)], shape=s))  
      })) # a `1 x n_test` matrix
      pred_y_draws[i, ] <- new_y
    }
  }
  else { # if mod="abs" 
    if (is.null(X_b_new)) X_b_new <- matrix(1, nrow=n_test, ncol=1) # Default design matrix for b
    if (ncol(X_b_new) != ncol(X_b)) stop("Dimensions of X_b_new and X_b must match.") 
    if (is.null(X_s_new)) X_s_new <- matrix(1, nrow=n_test, ncol=1) # Default design matrix for s
    if (ncol(X_s_new) != ncol(X_s)) stop("Dimensions of X_s_new and X_s must match.") 
    beta_b_ind <- grep("beta_b", colnames(parameter_draws)) # indices of betas for b
    beta_s_ind <- grep("beta_s", colnames(parameter_draws)) # indices of betas for s
    for (i in 1:n_draw){
      a <- parameter_draws[i, 1:n_train]
      logb <- parameter_draws[i, (n_train+1):(2*n_train)]
      beta_a <- parameter_draws[i, beta_a_ind]
      beta_b <- parameter_draws[i, beta_b_ind]
      beta_s <- parameter_draws[i, beta_s_ind]
      s <- s_draw_fun(i, (2*n_train+1):(3*n_train))
      X_all_a <- rbind(X_a, X_a_new)
      X_all_b <- rbind(X_b, X_b_new)
      X_all_s <- rbind(X_s, X_s_new)
      hyperparam_a1 <- exp(parameter_draws[i, "log_sigma_a"])
      hyperparam_b1 <- exp(parameter_draws[i, "log_sigma_b"])
      hyperparam_s1 <- exp(parameter_draws[i, "log_sigma_s"])
      # Construct conditional distribution function for a and logb
      if (kernel == "exp"){
	hyperparam_a2 <- exp(parameter_draws[i, "log_ell_a"])
	hyperparam_b2 <- exp(parameter_draws[i, "log_ell_b"])
	hyperparam_s2 <- exp(parameter_draws[i, "log_ell_s"])
	a_sim_fun <- sim_cond_normal(X_all_a%*%beta_a, a = a, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_exp, 
				     sigma = hyperparam_a1, 
				     ell = hyperparam_a2)
	logb_sim_fun <- sim_cond_normal(X_all_b%*%beta_b, a = logb, 
					locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
					kernel = kernel_exp, 
					sigma = hyperparam_b1, 
					ell = hyperparam_b2)
	s_sim_fun <- sim_cond_normal(X_all_s%*%beta_s, a = s, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_exp, 
				     sigma = hyperparam_s1, 
				     ell = hyperparam_s2)
      }
      else{
	hyperparam_a2 <- exp(parameter_draws[i, "log_kappa_a"])
	hyperparam_b2 <- exp(parameter_draws[i, "log_kappa_b"])
	hyperparam_s2 <- exp(parameter_draws[i, "log_kappa_s"])
        if (kernel == "spde"){
	  meshidxloc <- model$meshidxloc
          X_all_a <- rbind(as.matrix(X_a[meshidxloc,]), X_a_new)
          X_all_b <- rbind(as.matrix(X_b[meshidxloc,]), X_b_new)
          X_all_s <- rbind(as.matrix(X_s[meshidxloc,]), X_s_new)
	}
	a_sim_fun <- sim_cond_normal(X_all_a%*%beta_a, a = a, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_matern, sigma = hyperparam_a1, 
				     kappa = hyperparam_a2, nu = nu)
	logb_sim_fun <- sim_cond_normal(X_all_b%*%beta_b, a = logb, 
					locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
					kernel = kernel_matern, sigma = hyperparam_b1, 
					kappa = hyperparam_b2, nu = nu)
	s_sim_fun <- sim_cond_normal(X_all_s%*%beta_s, a = s, 
				     locs_new = as.matrix(locs_new), locs_obs = as.matrix(locs_obs), 
				     kernel = kernel_matern, 
				     sigma = hyperparam_s1, 
				     kappa = hyperparam_s2, nu = nu)
      }
      new_a <- a_sim_fun(1) # 1 x n_test matrix
      new_logb <- logb_sim_fun(1) # 1 x n_test matrix
      new_s <- s_sim_fun(1) # 1 x n_test matrix
      new_abs <- cbind(new_a, exp(new_logb), new_s) # A `1 x (3*n_test)` matrix
      new_y <- t(apply(X = new_abs, MARGIN = 1, FUN = function(row){
        unlist(Map(rgev, n=1, loc=row[1:n_test], scale=row[(n_test+1):(2*n_test)], 
		   shape=row[(2*n_test+1):length(row)]))  
      })) # a `1 x n_test` matrix
      pred_y_draws[i, ] <- new_y
    }
  }
  out <- list(pred_y_draws=pred_y_draws, locs_new=locs_new, locs_obs=locs_obs)
  class(out) <- "spatialGEVpred"
  out
}
