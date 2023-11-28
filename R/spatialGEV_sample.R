#' Get posterior parameter draws from a fitted GEV-GP model.
#'
#' @param model A fitted spatial GEV model object of class `spatialGEVfit`
#' @param n_draw Number of draws from the posterior distribution
#' @param observation whether to draw from the posterior distribution of the GEV observation?
#' @param loc_ind A vector of location indices to sample from. Default is all locations.
#' @return An object of class `spatialGEVsam`, which is a list of matrices containing the 
#' joint posterior draws of the parameters and optionally the GEV observations.
#' @examples 
#' \donttest{
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' n_loc <- nrow(locs)
#' beta_a <- mean(a); beta_b <- mean(logb)
#' fit <- spatialGEV_fit(y = y, locs = locs, random = "ab",
#'                       init_param = list(beta_a = beta_a,
#'                                         beta_b = beta_b,
#'                                         a = rep(0, n_loc), 
#'                                         log_b = rep(0, n_loc), 
#'                                         s = 0,
#'                                         log_sigma_a = 0, 
#'                                         log_kappa_a = 0,
#'                                         log_sigma_b = 0, 
#'                                         log_kappa_b = 0),
#'                       reparam_s = "positive",
#'                       kernel = "matern",
#'                       silent = TRUE) 
#' sam <- spatialGEV_sample(model = fit, n_draw = 100, 
#'                          observation = TRUE, loc_ind=1:10)
#' print(sam)
#' summary(sam)
#' }
#' @export
spatialGEV_sample <- function(model, n_draw, observation=FALSE, loc_ind=NULL){
  # Extract info from model
  kernel <- model$kernel
  rep <- model$report
  random <- model$random
  random_ind <- rep$env$random #indices of random effects
  n_loc <- length(unique(model$adfun$env$data$loc_ind)) # number of locations
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  ##########################
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
    else if (mod == "ab") {
      ind2rm <- setdiff(1:length(random_ind), c(meshidxloc, length(random_ind)/2 + meshidxloc))
      mean_random <- mean_random[-ind2rm]
      joint_cov <- joint_cov[-ind2rm, -ind2rm] 
    }
    else {# if mod="abs"
      ind2rm <- setdiff(1:length(random_ind), 
			c(meshidxloc, # indices of a in the random effects vec
			  length(random_ind)/3 + meshidxloc, # indices of log_b in the random vec
			  length(random_ind)/3*2 + meshidxloc)) # indices of s in  the random vec
      mean_random <- mean_random[-ind2rm]
      joint_cov <- joint_cov[-ind2rm, -ind2rm] 
    }
  }
  par_names_random <- paste0(names(mean_random), 1:n_loc) # add location index to the end of each parameter name
  par_names_fixed <- names(mean_fixed) # extract parameter names for the fixed effects
  joint_mean <- c(mean_random, mean_fixed)
  names(joint_mean) <- c(par_names_random, par_names_fixed) # modify parameter names
  fixed_ind <- (length(mean_random)+1):length(joint_mean) # indices of the fixed effects

  # Determine the positions of the samples based on location indices to be sampled
  if (is.null(loc_ind)){
    loc_ind <- 1:n_loc
    sam_ind <- 1:length(joint_mean)
  }
  else{
    loc_ind <- sort(loc_ind)
    patterns <- c(paste0("a", loc_ind), paste0("log_b", loc_ind), paste0("s", loc_ind))
    sam_ind <- which(names(joint_mean) %in% patterns)
    sam_ind <- c(sam_ind, fixed_ind) 
  }
  # Sample from joint normal
  joint_post_draw <- mvtnorm::rmvnorm(n_draw, joint_mean[sam_ind], joint_cov[sam_ind, sam_ind])
  output_list <- list(parameter_draws=joint_post_draw)
  # run below only if posterior draws of GEV observations are wanted
  if (observation){
    if (reparam_s == 3){ # unconstrained s
      s_fun <- function(x) x
    }
    else if (reparam_s == 1){ # positive s
      s_fun <- function(x) exp(x)
    }
    else if (reparam_s == 2){ # negative s
      s_fun <- function(x) -exp(x)
    }
    else{ # s=0
      s_fun <- function(x) 0
    }
    # Get parameter indices
    a_ind <- which(colnames(joint_post_draw)%in% paste0("a", loc_ind))
    if(mod=="a"){
      # b is fixed
      logb_ind <- "log_b"
    } else{
      # b is random
      logb_ind <- which(colnames(joint_post_draw)%in% paste0("log_b", loc_ind))
    }
    if(mod=="abs"){
      # s is random
      s_ind <- which(colnames(joint_post_draw)%in% paste0("s", loc_ind))
    } else{
      # s is fixed
      s_ind <- "s"
    }                   
    y_draws <- apply(joint_post_draw, 1,
	       function(all_draw){
		 mapply(rgev, n=1, loc=all_draw[a_ind], scale=exp(all_draw[logb_ind]),
		      shape=s_fun(all_draw[s_ind]))
	       })
    y_draws <- t(y_draws)
    colnames(y_draws) <- paste0("y", loc_ind)
    output_list[["y_draws"]] <- y_draws
  }
  class(output_list) <- "spatialGEVsam"
  output_list
}
