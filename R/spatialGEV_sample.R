#' Get posterior parameter draws from a fitted GEV-GP model.
#'
#' @param model A fitted spatial GEV model object of class `spatialGEVfit`
#' @param n_draw Number of draws from the posterior distribution
#' @param observation whether to draw from the posterior distribution of the GEV observation?
#' @param loc_ind A vector of location indices to sample from. Default is all locations.
#' @return An object of class `spatialGEVsam`, which is a list with the following elements:
#' \describe{
#'   \item{`parameter_draws`}{A matrix of joint posterior draws for the hyperparameters and the random effects at the `loc_ind` locations.}
#'   \item{`y_draws`}{If `observation == TRUE`, a matrix of corresponding draws from the posterior predictive GEV distribution at the `loc_ind` locations.}
#' }
#' @example examples/spatialGEV_sample.R
#' @export
spatialGEV_sample <- function(model, n_draw, observation=FALSE, loc_ind=NULL) {
  # Extract info from model
  rep <- model$report
  random <- model$random
  n_loc <- length(unique(model$adfun$env$data$loc_ind)) # number of locations
  reparam_s <- model$adfun$env$data$reparam_s # parametrization of s
  if(is.null(loc_ind)) loc_ind <- seq_len(n_loc)
  loc_ind <- seq_len(n_loc) %in% loc_ind # convert to logical
  # which parameters to keep in output
  sample_ind <- format_sample(adfun = model$adfun,
                              random = random,
                              loc_ind = loc_ind,
                              meshidxloc = model$meshidxloc)
  # construct mean vector
  mean_random <- rep$par.random
  mean_fixed <- rep$par.fixed
  mean_joint <- setNames(rep(NA, length(sample_ind)), names(sample_ind))
  mean_joint[names(mean_joint) %in% names(mean_random)] <- mean_random
  mean_joint[names(mean_joint) %in% names(mean_fixed)] <- mean_fixed
  # sampling
  prec_joint <- rep$jointPrecision
  if(!all(sapply(dimnames(prec_joint),
                 function(x) identical(x, names(mean_joint))))) {
    stop("Dimension name mismatch between `mean_joint` and `prec_joint`. Please file a bug report.")
  }
  joint_post_draw <- rmvn_prec(n_draw,
                               mean = mean_joint, prec = prec_joint)
  joint_post_draw <- joint_post_draw[,sample_ind]
  if(observation) {
    tmp_names <- names(sample_ind)[sample_ind]
    y_draw <- rgev_reparam(
      n = n_draw * sum(loc_ind),
      a = joint_post_draw[,tmp_names == "a"],
      log_b = joint_post_draw[,tmp_names == "log_b"],
      s = joint_post_draw[,tmp_names == "s"],
      reparam_s = reparam_s
    )
    y_draw <- matrix(y_draw, n_draw, sum(loc_ind))
  }
  # naming
  tmp_names <- colnames(joint_post_draw)
  for(random_nm in random) {
    tmp_names[tmp_names == random_nm] <- paste0(random_nm,
                                                which(loc_ind))
  }
  colnames(joint_post_draw) <- tmp_names
  output_list <- list(parameter_draws = joint_post_draw)
  if(observation) {
    colnames(y_draw) <- paste0("y", which(loc_ind))
    output_list$y_draws <- y_draw
  }
  class(output_list) <- "spatialGEVsam"
  output_list
}

#--- helper functions ----------------------------------------------------------

#' Get indices of random locations.
#'
#' @param adfun The `adfun` element of `model`.
#' @param random The `random` element of `model`.
#' @param loc_ind The location indices as a logical vector.
#' @param meshidxloc The `meshidxloc` element of the `model`, which can be `NULL`.
#' @return Named logical vector for which `parameter` elements should be sampled.  See Details.
#' @details The `parameter` elements -- be they fixed or random -- are returned in the order as they are specified by [TMB::MakeADFun()].  So the only thing that changes here is that some of the random effects are removed.
#' @noRd
format_sample <- function(adfun, random, loc_ind, meshidxloc) {
  sample_id <- rep(TRUE, length(adfun$env$par))
  sample_nm <- names(adfun$env$par)
  if(is.null(meshidxloc)) meshidxloc <- seq(1, length(loc_ind))
  for(random_nm in random) {
    # indices of given random effect
    random_id <- sample_nm == random_nm
    # determine which of these are excluded
    exclude_id <- rep(TRUE, sum(random_id))
    exclude_id[meshidxloc[loc_ind]] <- FALSE
    sample_id[which(random_id)[exclude_id]] <- FALSE
  }
  setNames(sample_id, sample_nm)
}

#' Sample from a multivariate normal with sparse precision matrix.
#'
#' @param n Number of random draws.
#' @param mean Mean vector.
#' @param prec Sparse precision matrix, i.e., inheriting from [Matrix::sparseMatrix] or its Cholesky factor, i.e., inheriting from [Matrix::CHMfactor-class].
#'
#' @return A matrix with `n` rows, each of which is a draw from the corresponding normal distribution.
#'
#' @details If the matrix is provided in precision form, it is converted to Cholesky form using `Matrix::Cholesky(prec, super = TRUE)`.  Once it is of form [Matrix::CHMfactor-class], this function is essentially copied from local function `rmvnorm()` in function `MC()` defined in [TMB::MakeADFun()].
#' @noRd
rmvn_prec <- function(n, mean, prec) {
  d <- ncol(prec) # number of mvn dimensions
  if(!is(prec, "CHMfactor")) {
    prec <- Matrix::Cholesky(prec, super = TRUE)
  }
  u <- matrix(rnorm(d*n),d,n)
  u <- Matrix::solve(prec,u,system="Lt")
  u <- Matrix::solve(prec,u,system="Pt")
  u <- t(as(u, "matrix") + mean)
}

#' Sample from GEV distribution with reparametrization.
#'
#' @param a Vector of location parameters.
#' @param log_b Vector of log-scale parameters.
#' @param s Vector of transformed shape parameters.  Ignored if `reparam_s == 0`.
#' @param reparam_s Integer type of shape parametrization.
#'
#' @details Largely copied from [evd::rgev()], except allows one to vectorize through the shape parameter, since `reparam_s` separately handles the special case `shape = 0`.  Thus, if `reparam_s =
#'
#' @noRd
rgev_reparam <- function(n, a, log_b, s, reparam_s) {
  loc <- a
  scale <- exp(log_b)
  if(reparam_s == 3) { # unconstrained s
    shape <- s
  } else if(reparam_s == 1) { # positive s
    shape <- exp(s)
  } else if(reparam_s == 2) { # negative s
    shape <- -exp(s)
  } else if(reparam_s == 0) { # s=0
    # sample from Gumbel distribution
    return(loc - scale * log(rexp(n)))
  } else {
    stop("Invalid value of `reparam_s`.")
  }
  return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}
