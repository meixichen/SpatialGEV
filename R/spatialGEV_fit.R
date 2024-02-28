#' Fit a GEV-GP model.
#'
#' @param data If `method == "laplace"`, a list of length `n_loc` where each
#'   element contains the GEV observations at the given spatial location.
#'   If `method == "maxsmooth"` as list with two elements: `est`,
#'   an `n_loc x 3` matrix of parameter estimates at each location,
#'   and `var`, a `3 x 3 x n_loc` array of corresponding variance estimates.
#' @param locs An `n_loc x 2` matrix of longitude and latitude of the corresponding response values.
#' @param random Either "a", "ab", or "abs", where `a` indicates the location parameter,
#' `b` indicates the scale parameter, `s` indicates the shape parameter.  This tells the model
#' which GEV parameters are considered as random effects.
#' @param method Either "laplace" or "maxsmooth".
#' @param init_param A list of initial parameters. See details.
#' @param reparam_s A flag indicating whether the shape parameter is "zero", "unconstrained",
#' constrained to be "negative", or constrained to be "positive". If model "abs" is used,
#' `reparam_s` cannot be zero. See details.
#' @param kernel Kernel function for spatial random effects covariance matrix. Can be "exp"
#' (exponential kernel), "matern" (Matern kernel), or "spde" (Matern kernel with SPDE
#' approximation described in Lindgren el al. 2011). To use the SPDE approximation,
#' the user must first install the INLA R package.
#' @param X_a `n_loc x r_a` design matrix for a, where `r-1` is the number of covariates. If not
#' provided, a `n_loc x 1` column matrix of 1s is used.
#' @param X_b `n_loc x r_b` design matrix for log(b). Does not need to be provided if b is fixed.
#' @param X_s `n_loc x r_s` design matrix for g(s), where g() is a transformation function of `s`.
#' Does not need to be provided if s is fixed.
#' @param nu Hyperparameter of the Matern kernel. Default is 1.
#' @param s_prior Optional. A length 2 vector where the first element is the mean of the normal
#' prior on s or log(s) and the second is the standard deviation. Default is NULL, meaning a
#' uniform prior is put on s if s is fixed, or a GP prior is applied if s is a random effect.
#' @param beta_prior Optional named list that specifies normal priors on the GP mean function
#' coefficients `beta`s. Each element of the list should be a named length 2 vector in which
#' the first element is mean and second element is sd.
#' E.g. `beta_prior=list(beta_a=c(0,100), beta_b=c(0,10), beta_s=c(-2,5))`.
#' Default is NULL, which means imposing a noninformative uniform flat prior.
#' @param matern_pc_prior Optional named list that specifies Penalized complexity
#' priors on the GP Matern covariance hyperparameters `sig` and `rho`, where `sig =
#' sqrt(sigma)` and `rho = sqrt(8*nu)/kappa`. Names must be `matern_a`, `matern_b`,
#' or `matern_s`.
#' E.g. `matern_pc_prior=list(matern_s=matern_pc_prior(100, 0.9, 2, 0.1))`.
#' Default is NULL, which means a flat prior. See `?matern_pc_prior` for more details.
#' @param return_levels Optional vector of return-level probabilities.
#' If provided, the posterior mean and standard deviation of the upper-tail GEV quantile at each
#' spatial location for each of these probabilities will be included in the summary output.
#' See `?summary.spatialGEV_fit` for details.
#' @param get_return_levels_cov Default is TRUE if `return_levels` is specified. Can be turned off
#' for when the number of locations is large so that the high-dimensional covariance matrix for
#' the return levels is not stored.
#' @param sp_thres Optional. Thresholding value to create sparse covariance matrix. Any distance
#' value greater than or equal to `sp_thres` will be set to 0. Default is -1, which means not
#' using sparse matrix. Caution: hard thresholding the covariance matrix often results in bad
#' convergence.
#' @param adfun_only Only output the ADfun constructed using TMB? If TRUE, model fitting is not
#' performed and only a TMB tamplate `adfun` is returned (along with the created mesh if kernel is
#' "spde").
#' This can be used when the user would like to use a different optimizer other than the default
#' `nlminb`. E.g., call `optim(adfun$par, adfun$fn, adfun$gr)` for optimization.
#' @param ignore_random Ignore random effect? If TRUE, spatial random effects are not integrated
#' out in the model. This can be helpful for checking the marginal likelihood.
#' @param silent Do not show tracing information?
#' @param mesh_extra_init A named list of scalars. Used when the SPDE kernel is used. The list
#' provides the initial values for a, log(b), and s on the extra triangles created in the mesh.
#' The default is `list(a=1, log_b=0, s=0.001)`.
#' @param get_hessian Default to TRUE so that `spatialGEV_sample()` can be used for sampling
#' from the Normal approximated posterior with the inverse Hessian as the Normal covariance.
#' @param ... Arguments to pass to `INLA::inla.mesh.2d()`. See details `?inla.mesh.2d()` and
#' Section 2.1 of Lindgren & Rue (2015) JSS paper.
#' This is used specifically for when `kernel="spde"`, in which case a mesh needs to be
#' constructed on the spatial domain. When no arguments are passed to `inla.mesh.2d()`, a
#' default argument is `max.edge=2`, which simply specifies the largest allowed triangle edge
#' length. It is strongly suggested that the user should specify these arguments if they would
#' like to use the SPDE kernel. Please make sure INLA package is installed before
#' using the SPDE approximation.
#' @return If `adfun_only=TRUE`, this function outputs a list returned by `TMB::MakeADFun()`.
#' This list contains components `par, fn, gr` and can be passed to an R optimizer.
#' If `adfun_only=FALSE`, this function outputs an object of class `spatialGEVfit`, a list
#; containing the following:
#' - An adfun object
#' - A fit object given by calling `nlminb()` on the adfun
#' - An object of class `sdreport` from TMB which contains the point estimates, standard error,
#' and precision matrix for the fixed and random effects
#' - Other helpful information about the model: kernel, data coordinates matrix, and optionally
#' the created mesh if `kernel="spde" (See details).
#'
#' `spatialGEV_model()` is used internally by `spatialGEV_fit()` to parse its inputs.  It returns a list with elements `data`, `parameters`, `random`, and `map` to be passed to [TMB::MakeADFun()].  If `kernel == "spde"`, the list also contains an element `mesh`.
#'
#' @details
#' This function adopts Laplace approximation using TMB model to integrate out the random effects.
#'
#' The random effects are assumed to follow Gaussian processes with mean 0 and covariance matrix
#' defined by the chosen kernel function. E.g., using the exponential kernel function:
#' ```
#' cov(i,j) = sigma*exp(-|x_i - x_j|/ell)
#' ```
#' When specifying the initial parameters to be passed to `init_param`, care must be taken to
#' count the number of parameters. Described below is how to specify `init_param` under different
#' settings of `random` and `kernel`. Note that the order of the parameters must match the
#' descriptions below (initial values specified below such as 0 and 1 are only examples).
#'
#' - random = "a", kernel = "exp":
#' `a` should be a vector and the rest are scalars. `log_sigma_a` and `log_ell_a` are
#' hyperparameters in the exponential kernel for the Gaussian process describing the spatial
#' variation of `a`.
#' ```
#' init_param = list(a = rep(1,n_locations), log_b = 0, s = 1,
#'                   beta_a = rep(0, n_covariates),
#'                   log_sigma_a = 0, log_ell_a = 0)
#' ```
#' Note that even if `reparam_s=="zero"`, an initial value for `s` still must be provided, even
#' though in this case the value does not matter anymore.
#'
#' - random = "ab", kernel = "exp":
#' When `b` is considered a random effect, its corresponding GP hyperparameters `log_sigma_b`
#' and `log_ell_b` need to be specified.
#' ```
#' init_param = list(a = rep(1,n_locations),
#'                   log_b = rep(0,n_locations), s=1,
#'                   beta_a = rep(0, n_covariates), beta_b = rep(0, n_covariates),
#'                   log_sigma_a = 0,log_ell_a = 0,
#'                   log_sigma_b = 0,log_ell_b = 0).
#' ```
#'
#' - random = "abs", kernel = "exp":
#' ```
#' init_param = list(a = rep(1,n_locations),
#'                   log_b = rep(0,n_locations),
#'                   s = rep(0,n_locations),
#'                   beta_a = rep(0, n_covariates),
#'                   beta_b = rep(0, n_covariates),
#'                   beta_s = rep(0, n_covariates),
#'                   log_sigma_a = 0,log_ell_a = 0,
#'                   log_sigma_b = 0,log_ell_b = 0).
#'                   log_sigma_s = 0,log_ell_s = 0).
#' ```
#'
#' - random = "abs", kernel = "matern" or "spde":
#' When the Matern or SPDE kernel is used, hyperparameters for the GP kernel are `log_sigma_a/b/s`
#' and `log_kappa_a/b/s` for each spatial random effect.
#' ```
#' init_param = list(a = rep(1,n_locations),
#'                   log_b = rep(0,n_locations),
#'                   s = rep(0,n_locations),
#'                   beta_a = rep(0, n_covariates),
#'                   beta_b = rep(0, n_covariates),
#'                   beta_s = rep(0, n_covariates),
#'                   log_sigma_a = 0,log_kappa_a = 0,
#'                   log_sigma_b = 0,log_kappa_b = 0).
#'                   log_sigma_s = 0,log_kappa_s = 0).
#' ```
#'
#' `raparam_s` allows the user to reparametrize the GEV shape parameter `s`. For example,
#' - if the data is believed to be right-skewed and lower bounded, this means `s>0` and one should
#' use `reparam_s = "positive"`;
#' - if the data is believed to be left-skewed and upper bounded, this means `s<0` and one should
#' use `reparam_s="negative"`.
#' - When `reparam_s = "zero"`, the data likelihood is a Gumbel distribution. In this case the data
#' has no upper nor lower bound. Finally, specify `reparam_s = "unconstrained"` if no sign
#' constraint should be imposed on `s`.
#'
#' Note that when reparam_s = "negative" or "postive", the initial value of `s` in `init_param`
#' should be that of log(|s|).
#'
#' When the SPDE kernel is used, a mesh on the spatial domain is created using
#' `INLA::inla.mesh.2d()`, which extends the spatial domain by adding additional triangles in the
#' mesh to avoid boundary effects in estimation. As a result, the number of `a` and `b`  will be
#' greater than the number of locations due to these additional triangles: each of them also has
#' their own `a` and `b` values. Therefore, the fit function will return a vector `meshidxloc` to
#' indicate the positions of the observed coordinates in the random effects vector.
#' @examples
#' \donttest{
#' library(SpatialGEV)
#' n_loc <- 20
#' a <- simulatedData$a[1:n_loc]
#' logb <- simulatedData$logb[1:n_loc]
#' logs <- simulatedData$logs[1:n_loc]
#' y <- simulatedData$y[1:n_loc]
#' locs <- simulatedData$locs[1:n_loc,]
#' # No covariates are included, only intercept is included.
#' fit <- spatialGEV_fit(y, locs = locs, random = "ab",
#'                       init_param = list(a = rep(0, n_loc),
#'                                         log_b = rep(0, n_loc),
#'                                         s = 0,
#'                                         beta_a = 0,
#'                                         beta_b = 0,
#'                                         log_sigma_a = 0,
#'                                         log_kappa_a = 0,
#'                                         log_sigma_b = 0,
#'                                         log_kappa_b = 0),
#'                       reparam_s = "positive",
#'                       kernel = "matern",
#'                       X_a = matrix(1, nrow=n_loc, ncol=1),
#'                       X_b = matrix(1, nrow=n_loc, ncol=1),
#'                       silent = TRUE)
#' print(fit)
#'
#' # To use a different optimizer other than the default `nlminb()`, create
#' # an object ready to be passed to optimizer functions using `adfun_only=TRUE`
#' obj <- spatialGEV_fit(y, locs = locs, random = "ab",
#'                       init_param = list(a = rep(0, n_loc),
#'                                         log_b = rep(0, n_loc),
#'                                         s = 0,
#'                                         beta_a = 0,
#'                                         beta_b = 0,
#'                                         log_sigma_a = 0,
#'                                         log_kappa_a = 0,
#'                                         log_sigma_b = 0,
#'                                         log_kappa_b = 0),
#'                       reparam_s = "positive",
#'                       kernel = "matern",
#'                       X_a = matrix(1, nrow=n_loc, ncol=1),
#'                       X_b = matrix(1, nrow=n_loc, ncol=1),
#'                       adfun_only = TRUE)
#' fit <- optim(obj$par, obj$fn, obj$gr)
#' }
#'
#' # Using the SPDE kernel (SPDE approximation to the Matern kernel)
#' # Make sure the INLA package is installed before using `kernel="spde"`
#' \dontrun{
#' library(INLA)
#' n_loc <- 20
#' y <- simulatedData2$y[1:n_loc]
#' locs <- simulatedData2$locs[1:n_loc,]
#' fit_spde <- spatialGEV_fit(y, locs = locs, random = "abs",
#'                            init_param = list(a = rep(0, n_loc),
#'                                              log_b = rep(0, n_loc),
#'                                              s = rep(-2, n_loc),
#'                                              beta_a = 0,
#'                                              beta_b = 0,
#'                                              beta_s = -2,
#'                                              log_sigma_a = 0,
#'                                              log_kappa_a = 0,
#'                                              log_sigma_b = 0,
#'                                              log_kappa_b = 0,
#'                                              log_sigma_s = 0,
#'                                              log_kappa_s = 0
#'                                              ),
#'                            reparam_s = "positive",
#'                            kernel = "spde",
#'                            beta_prior = list(beta_a=c(0,100), beta_b=c(0,10),
#'                                              beta_s=c(0,10)),
#'                            matern_pc_prior = list(
#'                                                  matern_a=matern_pc_prior(1e5,0.95,5,0.1),
#'                                                  matern_b=matern_pc_prior(1e5,0.95,3,0.1),
#'                                                  matern_s=matern_pc_prior(1e2,0.95,1,0.1)
#'                                                  ))
#' plot(fit_spde$mesh) # Plot the mesh
#' points(locs[,1], locs[,2], col="red", pch=16) # Plot the locations
#' }
#' @export
spatialGEV_fit <- function(data, locs, random = c("a", "ab", "abs"),
                           method = c("laplace", "maxsmooth"),
                           init_param, reparam_s,
                           kernel = c("spde", "matern", "exp"),
                           X_a = NULL, X_b = NULL, X_s = NULL, nu = 1,
                           s_prior = NULL, beta_prior = NULL,
                           matern_pc_prior = NULL,
                           return_levels=0., get_return_levels_cov=T,
                           sp_thres = -1, adfun_only = FALSE,
                           ignore_random = FALSE, silent = FALSE,
                           mesh_extra_init = list(a=0, log_b=-1, s=0.001),
                           get_hessian=TRUE,
                           ...) {
  # parse inputs
  kernel <- match.arg(kernel)
  random <- match.arg(random)
  method <- match.arg(method)
  if(method == "maxsmooth") {
    if((kernel != "spde") || (random != "abs")) {
      stop("For `method = 'maxsmooth'`, only `random = 'abs'` and `kernel = 'spde'` are currently implemented.")
    }
  }
  model <- spatialGEV_model(data = data, locs = locs, random = random,
                            method = method, init_param = init_param,
                            reparam_s = reparam_s, kernel = kernel,
                            X_a = X_a, X_b = X_b, X_s = X_s, nu = nu,
                            s_prior = s_prior, beta_prior = beta_prior,
                            matern_pc_prior = matern_pc_prior,
                            sp_thres = sp_thres, ignore_random = ignore_random,
                            mesh_extra_init = mesh_extra_init, ...)
  # Build TMB template
  model$data$return_levels <- return_levels
  adfun <- TMB::MakeADFun(data = model$data,
                          parameters = model$parameters,
                          random = model$random,
                          map = model$map,
                          DLL = "SpatialGEV_TMBExports",
                          silent = silent)
  # output
  if(adfun_only) {
    if(kernel == "spde") {
      out <- list(adfun = adfun, mesh = model$mesh)
    } else {
      out <- adfun
    }
  } else {
    start_t <- Sys.time()
    if (return_levels[1] != 0) {
      adfun_optim <- TMB::MakeADFun(data = c(model$data, return_levels=0.),
  				    parameters = model$parameters,
  				    random = model$random,
  				    map = model$map,
  				    DLL = "SpatialGEV_TMBExports",
  				    silent = silent)
      fit <- nlminb(adfun_optim$par, adfun_optim$fn, adfun_optim$gr)
      report <- TMB::sdreport(adfun, par.fixed = fit$par, 
                              getJointPrecision = get_hessian,
                              getReportCovariance = get_return_levels_cov)
    } else{
      adfun_optim <- adfun
      fit <- nlminb(adfun_optim$par, adfun_optim$fn, adfun_optim$gr)
      report <- TMB::sdreport(adfun_optim, getJointPrecision = get_hessian)
    }
    t_taken <- as.numeric(difftime(Sys.time(), start_t, units="secs"))
    out <- list(adfun = adfun_optim, fit = fit, report = report,
                time = t_taken, random = model$random, kernel = kernel,
                # FIXME: why not just call this locs?
                locs_obs = locs,
                X_a = model$data$design_mat_a,
                X_b = model$data$design_mat_b,
                X_s = model$data$design_mat_s,
                pdHess_avail = get_hessian & report$pdHess)
    if (return_levels[1] != 0) {
      n_probs <- length(return_levels)
      rl_inds <- lapply(1:n_probs, function(u) seq(u, length(report$value), by=n_probs))
      out_return_levels <- lapply(rl_inds, function(ind) report$value[ind])
      out_return_levels_sd <- lapply(rl_inds, function(ind) report$sd[ind])
      rl_list_names <- as.character(return_levels)
      names(out_return_levels) <- rl_list_names
      names(out_return_levels_sd) <- rl_list_names
      out$return_levels <- out_return_levels
      out$return_levels_sd <- out_return_levels_sd
      if (get_return_levels_cov){
        out_return_levels_cov <- lapply(rl_inds, function(ind) report$cov[ind, ind])
        names(out_return_levels_cov) <- rl_list_names
        out$return_levels_cov <- out_return_levels_cov
      }
    }
    if (kernel == "spde") {
      out$mesh <- model$mesh
      out$meshidxloc <- as.integer(model$mesh$idx$loc)
      out$nu <- nu
    } else if (kernel == "matern") {
      out$nu <- nu
    }
    class(out) <- "spatialGEVfit"
  }
  out
}
