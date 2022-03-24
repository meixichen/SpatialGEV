#' Fit a GEV-GP model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations. 
#' @param locs `n x 2` matrix of longitude and latitude of the corresponding response values.
#' @param random Either "a", "ab", or "abs", where `a` indicates the location parameter, 
#' `b` indicates the scale parameter, `s` indicates the shape parameter.  This tells the model
#' which GEV parameters are considered as random effects.
#' @param init_param A list of initial parameters. See details. 
#' @param reparam_s A flag indicating whether the shape parameter is "zero", "unconstrained", 
#' constrained to be "negative", or constrained to be "positive". If model "abs" is used, 
#' `reparam_s` cannot be zero. See details.
#' @param kernel Kernel function for spatial random effects covariance matrix. Can be "exp" 
#' (exponential kernel), "matern" (Matern kernel), or "spde" (Matern kernel with SPDE 
#' approximation described in Lindgren el al. 2011).
#' @param X_a `n x r` design matrix for a, where `r-1` is the number of covariates. If not 
#' provided, a `n x 1` column matrix of 1s is used.
#' @param X_b `n x r` design matrix for log(b). Does not need to be provided if b is fixed.
#' @param X_s `n x r` design matrix for g(s), where g() is a transformation function of `s`. 
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
#' @param ... Arguments to pass to `INLA::inla.mesh.2d()`. See details `?inla.mesh.2d()` and 
#' Section 2.1 of Lindgren & Rue (2015) JSS paper.
#' This is used specifically for when `kernel="spde"`, in which case a mesh needs to be 
#' constructed on the spatial domain. When no arguments are passed to `inla.mesh.2d()`, a 
#' default argument is `max.edge=2`, which simply specifies the largest allowed triangle edge
#' length. It is strongly suggested that the user should specify these arguments if they would 
#' like to use the SPDE kernel. 
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
#' \dontrun{
#' library(SpatialGEV)
#' a <- simulatedData$a
#' logb <- simulatedData$logb
#' logs <- simulatedData$logs
#' y <- simulatedData$y
#' locs <- simulatedData$locs
#' n_loc = nrow(locs)
#' # No covariates are included, only intercept is inlcuded.
#' fit <- spatialGEV_fit(y = y, locs = locs, random = "ab",
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
#' # Using the SPDE kernel (SPDE approximation to the Matern kernel)
#' fit_spde <- spatialGEV_fit(y = y, locs = locs, random = "abs",
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
#'                                              log_kappa_s = 0,
#'                                              ),
#'                            reparam_s = "positive",
#'                            kernel = "spde",
#'                            beta_prior = list(beta_a=c(0,100), beta_b=c(0,10),
#'                                              beta_s=c(0,10)),
#'                            matern_pc_prior = list(
#'                                                  matern_a=matern_pc_prior(1e5,0.95,5,0.1),
#'                                                  matern_b=matern_pc_prior(1e5,0.95,3,0.1),
#'                                                  matern_s=matern_pc_prior(1e2,0.95,1,0.1),
#'                                                  )
#'                            adfun_only = TRUE) 
#' library(INLA)
#' plot(fit_spde$mesh) # Plot the mesh
#' points(locs[,1], locs[,2], col="red", pch=16) # Plot the locations
#' }
#' @export
spatialGEV_fit <- function(y, locs, random, init_param, reparam_s, kernel="exp", 
			   X_a=NULL, X_b=NULL, X_s=NULL, nu=1, 
			   s_prior=NULL, beta_prior=NULL, matern_pc_prior=NULL, 
			   sp_thres=-1, adfun_only=FALSE, ignore_random=FALSE, 
			   silent=FALSE, 
			   mesh_extra_init=list(a=0, log_b=-1, s=0.001), ...){
  
  if (length(y) != nrow(locs)){
    stop("The length of y must be the same as the number of rows of locs.")
  }
  if (!(random %in% c("a", "ab", "abs"))){
    stop("Argument random must be either 'a', 'ab', or 'abs'.")
  } 
  if (reparam_s == "zero"){
    reparam_s <- as.integer(0)
    if (random == "abs") stop("When s is a random effect, reparam_s cannot be zero.")
  }
  else if (reparam_s == "positive"){
    reparam_s <- as.integer(1)
  }
  else if (reparam_s == "negative"){
    reparam_s <- as.integer(2)
  }
  else if (reparam_s == "unconstrained"){
    reparam_s <- as.integer(3)
  }
  else{
    stop("Argument reparam_s must be one of 'zero', 'unconstrained', 'positive', or 'negative'.")
  }  
  mod <- paste("model", random, sep="_")
  mod <- paste(mod, kernel, sep="_")
  n_loc <- length(y)
  n_obs <- sapply(y, length)
  if (missing(sp_thres)) sp_thres <- -1
  y <- unlist(y)

  #------ Prepare data input for TMB -------------
  if (kernel %in% c("exp", "matern")){
    # Default design matrices
    if (is.null(X_a)) X_a <- matrix(1, nrow=n_loc, ncol=1)
    if (is.null(X_b)) X_b <- matrix(1, nrow=n_loc, ncol=1)
    if (is.null(X_s)) X_s <- matrix(1, nrow=n_loc, ncol=1)
    
    dd <- as.matrix(stats::dist(locs))
    data <- list(model = mod, y = unlist(y), n_obs = n_obs, 
		 design_mat_a = X_a, design_mat_b = X_b, design_mat_s = X_s,
		 dd = dd, sp_thres = sp_thres, reparam_s = reparam_s)
    if (kernel == "matern") data$nu <- nu
  }
  else if (kernel == "spde"){
    mesh_args <- list(...)
    if (all(is.null(mesh_args$max.edge), 
	    is.null(mesh_args$max.n.strict), 
	    is.null(mesh_args$max.n))){ 
      # if none of the above is specified, use our default
      mesh <- INLA::inla.mesh.2d(locs, max.edge=2)
    } 
    else{
      mesh <- INLA::inla.mesh.2d(locs, ...)
    }
    spde <- (INLA::inla.spde2.matern(mesh)$param.inla)[c("M0", "M1", "M2")]
    n_s <- nrow(spde$M0) # number of mesh triangles created by INLA
    meshidxloc <- as.integer(mesh$idx$loc)
    
    # Default design matrices to provide to the data list
    if (is.null(X_a)){
      X_a <- matrix(1, nrow=n_s, ncol=1)
    }
    else{ # Expand the current design matrix using 0s due to the additional triangles in the mesh
      X_a_temp <- matrix(0, nrow=n_s, ncol=ncol(X_a))
      X_a_temp[,1] <- 1
      X_a_temp[meshidxloc,] <- X_a
      X_a <- X_a_temp   
    }
    if (is.null(X_b)){
      X_b <- matrix(1, nrow=n_s, ncol=1)
    }
    else{ # Expand the design matrix
      X_b_temp <- matrix(0, nrow=n_s, ncol=ncol(X_b))
      X_b_temp[,1] <- 1
      X_b_temp[meshidxloc,] <- X_b
      X_b <- X_b_temp   
    }
    # It is ok to have the additional element design_mat_b in the list even when it is not used 
    # in the TMB template
    data <- list(model = mod, y = unlist(y), n_obs = n_obs, 
	 	 design_mat_a = X_a, design_mat_b = X_b,  meshidxloc = meshidxloc-1, 
		 reparam_s = reparam_s, spde = spde, nu = nu)
    if (random == "a"){ 
      init_param_a <- rep(mesh_extra_init$a, n_s)
      init_param_a[meshidxloc] <- init_param$a # expand the vector of initial parameters due to extra location points introduced by mesh
      init_param$a <- init_param_a
    }
    else if (random == "ab") {
      init_param_a <- rep(mesh_extra_init$a, n_s)
      init_param_b <- rep(mesh_extra_init$log_b, n_s)
      init_param_a[meshidxloc] <- init_param$a 
      init_param_b[meshidxloc] <- init_param$log_b
      init_param$a <- init_param_a
      init_param$log_b <- init_param_b
    }
    else { # if random == "abs"
      if (is.null(X_s)){
	X_s <- matrix(1, nrow=n_s, ncol=1)
      }
      else{ # Expand the design matrix
	X_s_temp <- matrix(0, nrow=n_s, ncol=ncol(X_s))
	X_s_temp[,1] <- 1
	X_s_temp[meshidxloc,] <- X_s
	X_s <- X_s_temp   
      }
      data$design_mat_s <- X_s
      init_param_a <- rep(mesh_extra_init$a, n_s)
      init_param_b <- rep(mesh_extra_init$log_b, n_s)
      init_param_s <- rep(mesh_extra_init$s, n_s)
      init_param_a[meshidxloc] <- init_param$a 
      init_param_b[meshidxloc] <- init_param$log_b
      init_param_s[meshidxloc] <- init_param$s
      init_param$a <- init_param_a
      init_param$log_b <- init_param_b
      init_param$s <- init_param_s
    }
  }
  else {stop("kernel must be one of 'exp', 'matern', 'spde'!")}
  
  ############# Priors ##################### 
  # Optionally specify a normal prior on s if s is a fixed effect
  if (random != "abs"){ 
    if (missing(s_prior) | is.null(s_prior)){
      data$s_mean <- 9999
      data$s_sd <- 9999
    }
    else{ # specify the mean and sd of normal prior for s
      data$s_mean <- s_prior[1]
      data$s_sd <- s_prior[2]
    }
  }
  # Optionally specify normal priors on betas
  if (is.null(beta_prior)) {
    data$beta_prior <- as.integer(0)
    data$beta_a_prior <- rep(0,500)
    data$beta_b_prior <- rep(0,500)
    data$beta_s_prior <- rep(0,500)
  }
  else if (is.list(beta_prior)) {
    data$beta_prior <- as.integer(1)
    data$beta_a_prior <- beta_prior$beta_a
    data$beta_b_prior <- beta_prior$beta_b
    data$beta_s_prior <- beta_prior$beta_s 
  }
  else {
    stop("Check beta_prior.") 
  }
  # Optionally specify PC priors on Matern 
  if (kernel %in% c("matern", "spde")){
    data$a_pc_prior <- data$b_pc_prior <- data$s_pc_prior <- as.integer(0)
    data$range_a_prior <- data$range_b_prior <- data$range_s_prior <- c(1e5, 0.9)
    data$sigma_a_prior <- data$sigma_b_prior <- data$sigma_s_prior <- c(2, 0.1)
    if (!is.null(matern_pc_prior$matern_a)){
      if (class(matern_pc_prior$matern_a)=="PC_prior"){
	data$a_pc_prior <- as.integer(1)
	data$range_a_prior <- matern_pc_prior$matern_a$range_prior
	data$sigma_a_prior <- matern_pc_prior$matern_a$sigma_prior
      }
      else{
        stop("List elements must be provided using `matern_pc_prior()`")
      }
    }
    if (!is.null(matern_pc_prior$matern_b)){
      if (class(matern_pc_prior$matern_b)=="PC_prior"){
	data$b_pc_prior <- as.integer(1)
	data$range_b_prior <- matern_pc_prior$matern_b$range_prior
	data$sigma_b_prior <- matern_pc_prior$matern_b$sigma_prior
      }
      else{
        stop("List elements must be provided using `matern_pc_prior()`")
      }
    }
    if (!is.null(matern_pc_prior$matern_s)){
      if (class(matern_pc_prior$matern_s)=="PC_prior"){
	data$s_pc_prior <- as.integer(1)
	data$range_s_prior <- matern_pc_prior$matern_s$range_prior
	data$sigma_s_prior <- matern_pc_prior$matern_s$sigma_prior
      }
      else{
        stop("List elements must be provided using `matern_pc_prior()`")
      }
    }
    else if (!is.null(matern_pc_prior) & !is.list(matern_pc_prior)) {
      stop("Check matern_pc_prior: must be a named list with names one or more of
	   matern_a, matern_b, or matern_s, and the elements must be provided using
	   `matern_pc_prior()` function.")
    }
  }
  #------ End: prepare data input for TMB ----------------
  
  if (ignore_random){
    random <- NULL
  }
  else if (random == "ab" & !ignore_random){ 
    random <- c("a", "log_b")
  }
  else if (random == "abs" & !ignore_random){
    random <- c("a", "log_b", "s")
  } 
  #else random <- "a"

  # If using Gumbel, make sure s is not being estimated
  map <- list()
  if (reparam_s == "zero") { 
    map <- list(s = factor(NA))
  }

  # Build TMB template
  adfun <- TMB::MakeADFun(data = data,
                          parameters = init_param,
                          random = random,
                          map = map,
                          DLL = "SpatialGEV_TMBExports", 
                          silent = silent)
  
  if (adfun_only){
    if (kernel == "spde"){
      list(adfun=adfun, mesh=mesh)
    }
    else{
      adfun
    }
  }
  else{
    start_t <- Sys.time()
    fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
    report <- TMB::sdreport(adfun, getJointPrecision = TRUE)
    t_taken <- as.numeric(difftime(Sys.time(), start_t, unit="secs"))
    out <- list(adfun=adfun, fit=fit, report=report, 
		time=t_taken, random=random, kernel=kernel, 
		locs_obs=locs, X_a=X_a, X_b=X_b, X_s=X_s)
    if (kernel == "spde"){
      out$mesh <- mesh
      out$meshidxloc <- meshidxloc
      out$nu <- nu
    }
    else if (kernel == "matern"){
      out$nu <- nu
    }
    class(out) <- "spatialGEVfit"
    out 
  }
}
