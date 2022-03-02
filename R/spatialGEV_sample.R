#' Get posterior parameter draws from a fitted GEV-GP model.
#'
#' @param model A fitted spatial GEV model object of class `spatialGEVfit`
#' @param n_draw Number of draws from the posterior distribution
#' @param observation whether to draw from the posterior distribution of the GEV observation?
#' @return An object of class `spatialGEVsam`, which is a list of matrices containing the 
#' joint posterior draws of the parameters and optionally the GEV observations.
#' @examples 
#' \dontrun{
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' n_loc = nrow(locs)
#' fit <- spatialGEV_fit(y = y, X = locs, random = "ab",
#'                       init_param = list(a = rep(0, n_loc), 
#'                                         log_b = rep(0, n_loc), 
#'                                         s = 0,
#'                                         log_sigma_a = 0, 
#'                                         log_kappa_a = 0,
#'                                         log_sigma_b = 0, 
#'                                         log_kappa_b = 0),
#'                       reparam_s = "positive",
#'                       kernel = "matern",
#'                       silent = TRUE) 
#' sam <- spatialGEV_sample(model = fit, n_draw = 5000, observation = TRUE)
#' print(sam)
#' summary(sam)
#' }
#' @export
spatialGEV_sample <- function(model, n_draw, observation=FALSE){
  # Extract info from model
  kernel <- model$kernel
  rep <- model$report
  random_ind <- rep$env$random #indices of random effects
  n_loc <- length(model$adfun$env$data$n_obs) # number of locations
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  ##########################
  if (length(random_ind) == n_loc) {
    mod <- "a"
  }
  else if (length(random_ind) >= 2*n_loc){
    mod <- "ab"
  }
  else {
    stop("Cannot identify which GEV parameters are random.")
  }
  
  # Sample from MVN
  jointPrec_mat <- rep$jointPrecision
  C <- chol(jointPrec_mat)
  joint_cov <- backsolve(r = C, 
                         x = backsolve(r = C, x = diag(nrow(jointPrec_mat)), 
                                       transpose = TRUE)) # Cholesky decomp
  mean_random <- rep$par.random
  mean_fixed <- rep$par.fixed
  # Extract locations of the data if using SPDE (b/c there are additional boundary locations
  # included for fitting purpose
  if (kernel == "spde"){ 
    meshidxloc <- model$meshidxloc
    if (mod == "a") { 
      mean_random <- mean_random[meshidxloc]
      ind2rm <- setdiff(1:length(random_ind), meshidxloc) 
      joint_cov <- joint_cov[-ind2rm, -ind2rm]
    }
    else {
      ind2rm <- setdiff(1:length(random_ind), c(meshidxloc, length(random_ind)/2 + meshidxloc))
      mean_random <- mean_random[-ind2rm]
      joint_cov <- joint_cov[-ind2rm, -ind2rm] 
    }
  }
  par_names_random <- paste0(names(mean_random), 1:n_loc) # add location index to the end of each parameter name
  par_names_fixed <- names(mean_fixed) # extract parameter names for the fixed effects
  joint_mean <- c(mean_random, mean_fixed)
  names(joint_mean) <- c(par_names_random, par_names_fixed) # modify parameter names
  joint_post_draw <- mvtnorm::rmvnorm(n_draw, joint_mean, joint_cov) # sample from the joint MVN
  output_list <- list(parameter_draws=joint_post_draw)
  # run below only if posterior draws of GEV observations are wanted
  if (observation){
    s_draw_fun <- function(j, s_ind){ 
      # `j` points to the j-th draw and `s_ind` is the index of s in the parameter vector of the model
      if (reparam_s == 3){
        joint_post_draw[j, s_ind]
      }
      else if (reparam_s == 1){
        exp(joint_post_draw[j, s_ind])
      }
      else if (reparam_s == 2){
        -exp(joint_post_draw[j, s_ind])
      }
      else{ # if reparam_s=0
        0
      }
    }
    y_draws <- rep(NULL, n_draw) # <<==== TODO: Should allocate space for y_draws. Tried to preallocate a large matrix but got error message about not enough memory.
    for (i in 1:n_loc){
      loc_draw <- rep(NA, n_draw)
      for (j in 1:n_draw){
        a_draw <- joint_post_draw[j, i]
        if (mod == "a"){
          b_draw <- joint_post_draw[j, n_loc+1]
          s_draw <- s_draw_fun(j, n_loc+2)
        }
        else {
          b_draw <- joint_post_draw[j, i+n_loc]
          s_draw <- s_draw_fun(j, 2*n_loc+1)
        }
        loc_draw[j] <- rgev(1, loc = a_draw, scale = exp(b_draw), shape = s_draw)
      }
      y_draws <- cbind(y_draws, loc_draw)
    }
    colnames(y_draws) <- paste0("y", 1:n_loc)
    output_list[["y_draws"]] <- y_draws
  }
  class(output_list) <- "spatialGEVsam"
  output_list
}
