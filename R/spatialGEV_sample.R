#' Get posterior draws from a fitted spatial GEV model.
#'
#' @param model A fitted spatial GEV model returned by `spatialGEV_fit`
#' @param n_draw Number of draws from the posterior distribution
#' @param observation whether to draw from the posterior distribution of the GEV observation?
#' @return A list of matrices containing the joint posterior draws of the parameters and optionally the GEV observations
#' @export
spatialGEV_sample <- function(model, n_draw, observation=TRUE){
  # Extract info from model
  rep <- model$report
  random_ind <- rep$env$random #indices of random effects
  n_obs <- nrow(model$adfun$env$data$dd) # number of locations
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  ##########################
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
  joint_cov <- backsolve(r = C, 
                         x = backsolve(r = C, x = diag(nrow(jointPrec_mat)), 
                                       transpose = TRUE)) # Cholesky decomp
  mean_random <- rep$par.random
  mean_fixed <- rep$par.fixed
  par_names_random <- paste0(names(mean_random), 1:n_obs) # add location index to the end of each parameter name
  par_names_fixed <- names(mean_fixed) # extract parameter names for the fixed effects
  joint_mean <- c(mean_random, mean_fixed)
  names(joint_mean) <- c(par_names_random, par_names_fixed) # modify parameter names
  joint_post_draw <- mvtnorm::rmvnorm(n_draw, joint_mean, joint_cov) # sample from the joint MVN
  output_list <- list(parameter_draws=joint_post_draw)
  # run below only if posterior draws of GEV observations are wanted
  if (observation){
    s_draw_fun <- function(j, s_ind){ 
      # `j` points to the j-th draw and `s_ind` is the index of s in the parameter vector of the model
      if (reparam_s == "unconstrained"){
        joint_post_draw[j, s_ind]
      }
      else if (reparam_s == "positive"){
        exp(joint_post_draw[j, s_ind])
      }
      else if (reparam_s == "negative"){
        -exp(joint_post_draw[j, s_ind])
      }
      else{ # if reparam_s="zero"
        0
      }
    }
    y_draws <- rep(NULL, n_draw) # <<==== TODO: Should allocate space for y_draws. Tried to preallocate a large matrix but got error message about not enough memory.
    for (i in 1:n_obs){
      loc_draw <- rep(NA, n_draw)
      for (j in 1:n_draw){
        a_draw <- joint_post_draw[j, i]
        if (mod == "a"){
          b_draw <- joint_post_draw[j, n_obs+1]
          s_draw <- s_draw_fun(j, n_obs+2)
        }
        else {
          b_draw <- joint_post_draw[j, i+n_obs]
          s_draw <- s_draw_fun(j, 2*n_obs+1)
        }
        loc_draw[j] <- rgev(1, loc = a_draw, scale = exp(b_draw), shape = s_draw)
      }
      y_draws <- cbind(y_draws, loc_draw)
    }
    colnames(y_draws) <- paste0("y", 1:n_obs)
    output_list[["y_draws"]] <- y_draws
  }
  output_list
}
