#' @rdname spatialGEV_fit
#' @export
spatialGEV_model <- function(data, locs, random = c("a", "ab", "abs"),
                             method = c("laplace", "maxsmooth"),
                             init_param, reparam_s,
                             kernel = c("spde", "matern", "exp"),
                             X_a = NULL, X_b = NULL, X_s = NULL, nu = 1,
                             s_prior = NULL, beta_prior = NULL,
                             matern_pc_prior = NULL,
                             sp_thres = -1, ignore_random = FALSE,
                             mesh_extra_init = list(a=0, log_b=-1, s=0.001),
                             ...) {
  method <- match.arg(method)
  random <- parse_random(random)
  out_data <- parse_data(data, locs = locs, random = random, method = method)
  ## y <- data_out$y
  ## loc_ind <- data_out$loc_ind
  reparam_s <- parse_reparam_s(reparam_s, random = random)
  kernel <- match.arg(kernel)
  #------ Prepare data input for TMB -------------
  data <- list(model = parse_model(random = random,
                                   kernel = kernel, method = method),
               reparam_s = reparam_s)
  if(method == "laplace") {
    data <- c(data, list(y = out_data$y))
  } else if(method == "maxsmooth") {
    data <- c(data,
              list(obs = out_data$random_est,
                   cov_obs = out_data$random_var,
                   reparam_s = reparam_s))
  }
  if(kernel %in% c("exp", "matern")) {
    out_kernel <- parse_kernel_basic(locs = locs,
                                     X_a = X_a, X_b = X_b, X_s = X_s)
    data <- c(data,
              list(loc_ind = out_data$loc_ind-1,
                   design_mat_a = out_kernel$X_a,
                   design_mat_b = out_kernel$X_b,
                   design_mat_s = out_kernel$X_s,
                   dd = out_kernel$dd,
                   sp_thres = sp_thres))
    if(kernel == "matern") data$nu <- nu
  } else if(kernel == "spde") {
    out_kernel <- parse_kernel_spde(locs = locs, X_a = X_a, X_b = X_b, X_s,
                                    loc_ind = out_data$loc_ind,
                                    init_param = init_param, random = random,
                                    mesh_extra_init = mesh_extra_init, ...)
    # It is ok to have the additional element design_mat_b in the list
    # even when it is not used in the TMB template
    data <- c(data,
              list(loc_ind = out_kernel$loc_ind-1,
                   design_mat_a = out_kernel$X_a,
                   design_mat_b = out_kernel$X_b,
                   design_mat_s = out_kernel$X_s,
                   spde = out_kernel$spde,
                   nu = nu))
    init_param <- out_kernel$init_param
  }
  ############# Priors #####################
  out_priors <- parse_priors(random = random, kernel = kernel,
                             s_prior = s_prior,
                             beta_prior = beta_prior,
                             matern_pc_prior = matern_pc_prior)
  data <- c(data, out_priors)
  #------ End: prepare data input for TMB ----------------
  if(ignore_random) {
    random <- NULL
  } else {
    random <- names(random)[random]
  }
  # If using Gumbel, make sure s is not being estimated
  map <- list()
  if (reparam_s == 0L) {
    map <- list(s = factor(NA))
  }
  out <- list(data = data, parameters = init_param, random = random, map = map)
  if(kernel == "spde") out$mesh <- out_kernel$mesh
  out
}

#' @noRd
#'
#' @return For `method == "laplace"`, a list with elements `y`, `loc_ind`.  For `method == "maxsmooth"`, a list with elements `random_est`, `random_var`, `loc_ind`.
parse_data <- function(data, locs, random,
                       method = c("laplace", "maxsmooth")) {
  match.arg(method)
  n_loc <- nrow(locs)
  n_par <- sum(random)
  if(method == "laplace") {
    y <- data
    if(!is.list(y) ||
       !all(sapply(y, is.numeric))) {
      stop("For `method == 'laplace', `data must be a numeric list.")
    } else if(length(y) != n_loc) {
      stop("For `method == 'laplace', must have `length(data) == nrow(locs)`.")
    }
    ## n_loc <- length(y)
    n_obs <- sapply(y, length)
    y <- unlist(y)
    loc_ind <- rep(1:n_loc, times=n_obs) # location ind associated with each obs
    out <- list(y = y, loc_ind = loc_ind)
  } else if(method == "maxsmooth") {
    if(!all(c("est", "var") %in% names(data)) ||
       !is.numeric(data$est) ||
       !is.numeric(data$var)) {
      stop("For `method == 'maxsmooth'`, `data` must be a list with named elements 'est' and 'var'.")
    } else if(!isTRUE(all(dim(data$est) == c(n_loc, n_par)))) {
      stop("Incorrect dimensions for `data$est`.")
    } else if(!isTRUE(all(dim(data$var) == c(n_par, n_par, n_loc)))) {
      stop("Incorrect dimensions for `data$var`.")
    }
    out <- list(random_est = t(data$est),
                random_var = matrix(data$var, n_par, n_par * n_loc),
                loc_ind = 1:n_loc)
  }
  out
}

#' @noRd
#' @return A logical with names `a`, `log_b`, `s`.
parse_random <- function(random = c("a", "ab", "abs")) {
  random <- match.arg(random)
  if (!(random %in% c("a", "ab", "abs"))) {
    stop("Argument random must be either 'a', 'ab', or 'abs'.")
  }
  c(a = random %in% c("a", "ab", "abs"),
    log_b = random %in% c("ab", "abs"),
    s = random %in% "abs")
}

#' @noRd
#' @details FIXME: Better naming convention for TMB files.
parse_model <- function(random, kernel, method) {
  param_names <- paste0(c("a", "b", "s")[random], collapse = "")
  model <- paste("model",
                 param_names,
                 kernel,
                 sep = "_")
  if(method == "maxsmooth") model <- paste0(model, "_maxsmooth")
  model
}


#' @noRd
#' @return An integer between 1 and 4.
parse_reparam_s <- function(reparam_s, random) {
  if(reparam_s == "zero") {
    reparam_s <- as.integer(0)
    if (random["s"]) stop("When s is a random effect, reparam_s cannot be zero.")
  } else if (reparam_s == "positive") {
    reparam_s <- as.integer(1)
  } else if (reparam_s == "negative") {
    reparam_s <- as.integer(2)
  } else if (reparam_s == "unconstrained") {
    reparam_s <- as.integer(3)
  } else {
    stop("Argument reparam_s must be one of 'zero', 'unconstrained', 'positive', or 'negative'.")
  }
  reparam_s
}

#' @noRd
#' @return For `kernel %in% c("exp", "matern")`. A list with elements `X_a`, `X_b`, `X_s`, `dd`.
parse_kernel_basic <- function(locs, X_a, X_b, X_s) {
  n_loc <- nrow(locs)
  # Default design matrices
  if(is.null(X_a)) X_a <- matrix(1, nrow=n_loc, ncol=1)
  if(is.null(X_b)) X_b <- matrix(1, nrow=n_loc, ncol=1)
  if(is.null(X_s)) X_s <- matrix(1, nrow=n_loc, ncol=1)
  dd <- as.matrix(stats::dist(locs))
  out <- list(X_a = X_a, X_b = X_b, X_s = X_s, dd = dd)
  out
}

#' @noRd
#' @return A list with elements `X_a`, `X_b`, `X_s`, `spde`.
parse_kernel_spde <- function(locs, X_a, X_b, X_s,
                              loc_ind,
                              init_param, random, mesh_extra_init, ...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Please install package 'INLA' if using 'kernel='spde''.")
  }
  mesh_args <- list(...)
  if(all(is.null(mesh_args$max.edge),
         is.null(mesh_args$max.n.strict),
         is.null(mesh_args$max.n))) {
    # if none of the above is specified, use our default
    mesh <- INLA::inla.mesh.2d(locs, max.edge=2)
  } else {
    mesh <- INLA::inla.mesh.2d(locs, ...)
  }
  spde <- (INLA::inla.spde2.matern(mesh)$param.inla)[c("M0", "M1", "M2")]
  n_s <- nrow(spde$M0) # number of mesh triangles created by INLA
  meshidxloc <- as.integer(mesh$idx$loc)
  out <- lapply(list(X_a = X_a, X_b = X_b, X_s = X_s), function(X) {
    if (is.null(X)) {
      # Default design matrices to provide to the data list
      X <- matrix(1, nrow=n_s, ncol=1)
    } else {
      # Expand the current design matrix using 0s due to
      # the additional triangles in the mesh
      X_temp <- matrix(0, nrow=n_s, ncol=ncol(X))
      X_temp[,1] <- 1
      X_temp[meshidxloc,] <- X
      X <- X_temp
    }
    X
  })
  out$spde <- spde
  out$mesh <- mesh
  # expand init_param due to extra location points introduced by mesh
  for(nm in names(random)[random]) {
    param_new <- rep(mesh_extra_init[[nm]], n_s)
    param_new[meshidxloc] <- init_param[[nm]]
    init_param[[nm]] <- param_new
  }
  out$init_param <- init_param
  out$loc_ind <- meshidxloc[loc_ind]
  out
}

#' @noRd
#' @return A list with all prior elements.
parse_priors <- function(random, kernel,
                         s_prior, beta_prior, matern_pc_prior) {
  prior_list <- list()
  if(!random["s"]) {
    # specify a normal prior on s if s is a fixed effect
    if(is.null(s_prior)) {
      s_prior <- c(9999, 9999)
    }
    prior_list$s_mean <- s_prior[1]
    prior_list$s_sd <- s_prior[2]
  }
  # Optionally specify normal priors on betas
  if(is.null(beta_prior)) {
    prior_list$beta_prior <- as.integer(0)
    prior_list$beta_a_prior <- c(0,500)
    prior_list$beta_b_prior <- c(0,500)
    prior_list$beta_s_prior <- c(0,500)
  } else if (is.list(beta_prior)) {
    prior_list$beta_prior <- as.integer(1)
    prior_list$beta_a_prior <- beta_prior$beta_a
    prior_list$beta_b_prior <- beta_prior$beta_b
    prior_list$beta_s_prior <- beta_prior$beta_s
  } else {
    stop("Check beta_prior.")
  }
  # Optionally specify PC priors on Matern
  if(kernel %in% c("matern", "spde")) {
    if(!is.null(matern_pc_prior) && !is.list(matern_pc_prior)) {
      stop("Check matern_pc_prior: must be a named list with names one or more of
	   `matern_a`, `matern_b`, or `matern_s`, and the elements must be provided using the
	   `matern_pc_prior()` function.")
    }
    for(par in c("a", "b", "s")) {
      # names of list elements wrt par
      par_pc_prior <- paste0(par, "_pc_prior")
      range_par_prior <- paste0("range_", par, "_prior")
      sigma_par_prior <- paste0("sigma_", par, "_prior")
      matern_par <- paste0("matern_", par)
      # create prior
      prior_list[[par_pc_prior]] <- as.integer(0)
      prior_list[[range_par_prior]] <- c(1e5, 0.9)
      prior_list[[sigma_par_prior]] <- c(2, 0.1)
      if(!is.null(matern_pc_prior[[matern_par]])) {
        if(inherits(matern_pc_prior[[matern_par]], "PC_prior")) {
          prior_list[[par_pc_prior]] <- as.integer(1)
          prior_list[[range_par_prior]] <- matern_pc_prior[[matern_par]]$range_prior
          prior_list[[sigma_par_prior]] <- matern_pc_prior[[matern_par]]$sigma_prior
        } else {
          stop(paste0("`matern_pc_prior$", matern_par,
                      "` must be created with `matern_pc_prior()`."))
        }
      }
    }
  }
  prior_list
}
