#' Fit a GEV-GP model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations. 
#' @param X `n x 2` matrix of longitude and latitude of the corresponding response values.
#' @param random Either "a" or "ab". This indicates which GEV parameters are considered as 
#' random effects.
#' @param init_param A list of initial parameters. See details. 
#' @param reparam_s A flag indicating whether the shape parameter is "zero", "unconstrained", 
#' constrained to be "negative", or constrained to be "positive". See details.
#' @param kernel Kernel function for spatial random effects covariance matrix. Can be "exp" 
#' (exponential kernel), "matern" (Matern kernel), or "spde" (Matern kernel with SPDE 
#' approximation described in Lindgren el al. 2011).
#' @param nu Hyperparameter of the Matern kernel. Default is 1. 
#' @param s_prior Optional. A length 2 vector where the first element is the mean of the normal 
#' prior on s or log(s) and the second is the standard deviation. Default is NULL, meaning a 
#' uniform prior is put on s.
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
#' @param ... Arguments to pass to `INLA::inla.mesh.2d()`. See details `?inla.mesh.2d()` and 
#' Section 2.1 of Lindgren & Rue (2015) JSS paper.
#' This is used specifically for when `kernel="spde"`, in which case a mesh needs to be 
#' constructed on the spatial domain. When no arguments are passed to `inla.mesh.2d()`, a 
#' default argument is `max.edge=c(1,2)`, which simply specifies the largest allowed triangle edge
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
#' descriptions below.
#'
#' - random = "a", kernel = "exp": 
#' `a` should be a vector and the rest are scalars. `log_sigma_a` and `log_ell_a` are 
#' hyperparameters in the exponential kernel for the Gaussian process describing the spatial 
#' variation of `a`.  
#' ```
#' init_param = list(a = rep(1,n_locations), log_b = 0, s = 1,
#'                   log_sigma_a = 0, log_ell_a = 0)
#' ```
#'
#' - random = "ab", kernel = "exp": If
#' When the scale parameter `b` is considered a random effect, its corresponding GP hyperparameters 
#' `log_sigma_b` and `log_ell_b` need to be specified.
#' ```
#' init_param = list(a = rep(1,n_locations),
#'                   log_b = rep(0,n_locations), s=1,
#'                   log_sigma_a = 0,log_ell_a = 0, 
#'                   log_sigma_b = 0,log_ell_b = 0).
#' ```
#' 
#' - random = "ab", kernel = "matern" or "spde": 
#' When the Matern or SPDE kernel is used, hyperparameters for the GP kernel are `log_sigma_a/b` 
#' and `log_kappa_a/b` for each spatial random effect. 
#' ``` 
#' init_param = list(a = rep(1,n_locations),
#'                   log_b = rep(0,n_locations), s=1,
#'                   log_sigma_a = 0,log_kappa_a = 0, 
#'                   log_sigma_b = 0,log_kappa_b = 0).
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
#' `INLA::inla.mesh.2d()', which extends the spatial domain by adding additional triangles in the
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
#' print(fit)
#' 
#' # Using the SPDE kernel (SPDE approximation to the Matern kernel)
#' fit_spde <- spatialGEV_fit(y = y, X = locs, random = "ab",
#'                            init_param = list(a = rep(0, n_loc), 
#'                                              log_b = rep(0, n_loc), 
#'                                              s = 0,
#'                                              log_sigma_a = 0, 
#'                                              log_kappa_a = 0,
#'                                              log_sigma_b = 0, 
#'                                              log_kappa_b = 0),
#'                            reparam_s = "positive",
#'                            kernel = "spde",
#'                            adfun_only = TRUE) 
#' library(INLA)
#' plot(fit_spde$mesh) # Plot the mesh
#' points(locs[,1], locs[,2], col="red", pch=16) # Plot the locations
#' }
#' @export
spatialGEV_fit <- function(y, X, random, init_param, reparam_s, kernel="exp", nu=1, s_prior= NULL, sp_thres=-1, adfun_only=FALSE, ignore_random=FALSE, silent=FALSE, ...){
  
  if (length(y) != nrow(X)){
    stop("The length of y must be the same as the number of rows of X.")
  }
  if (!(random %in% c("a", "ab"))){
    stop("Argument random must be either 'a' or 'ab'.")
  } 
  if (reparam_s == "zero"){
    reparam_s <- as.integer(0)
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
    dd <- as.matrix(stats::dist(X))
    data <- list(model = mod, y = unlist(y), n_obs = n_obs, dd = dd, sp_thres = sp_thres, reparam_s = reparam_s)
    if (kernel == "matern") data$nu <- nu
  }
  else if (kernel == "spde"){
    mesh_args <- list(...)
    if (all(is.null(mesh_args$max.edge), 
	    is.null(mesh_args$max.n.strict), 
	    is.null(mesh_args$max.n))){ 
      # if none of the above is specified, use our default
      mesh <- INLA::inla.mesh.2d(X, max.edge=c(1,2))
    } 
    else{
      mesh <- INLA::inla.mesh.2d(X, ...)
    }
    spde <- (INLA::inla.spde2.matern(mesh)$param.inla)[c("M0", "M1", "M2")]
    n_s <- nrow(spde$M0) # number of mesh triangles created by INLA
    meshidxloc <- as.integer(mesh$idx$loc - 1)
    data <- list(model = mod, y = unlist(y), n_obs = n_obs, meshidxloc = meshidxloc, 
		 reparam_s = reparam_s, spde = spde, nu = nu)
    if (random == "a"){ 
      init_param_a <- rep(0, n_s)
      init_param_a[meshidxloc+1] <- init_param$a # expand the vector of initial parameters due to extra location points introduced by mesh
      init_param$a <- init_param_a
    }
    else {
      init_param_a <- rep(0, n_s)
      init_param_b <- rep(-1, n_s)
      init_param_a[meshidxloc+1] <- init_param$a 
      init_param_b[meshidxloc+1] <- init_param$log_b
      init_param$a <- init_param_a
      init_param$log_b <- init_param_b
    }
  }
  else {stop("kernel must be one of 'exp', 'matern', 'spde'!")}
  #------ End: prepare data input for TMB ----------------

  if (missing(s_prior) | is.null(s_prior)){
    data$s_mean <- 9999
    data$s_sd <- 9999
  }
  else{ # specify the mean and sd of normal prior for s
    data$s_mean <- s_prior[1]
    data$s_sd <- s_prior[2]
  }
  
  if (random == "ab" & !ignore_random){
    random <- c("a", "log_b")
  }
  else if (ignore_random){
    random <- NULL
  }
  
  map <- list()
  if (reparam_s == "zero") { # if using Gumbel, make sure s is not being estimated
    map <- list(s = factor(NA))
  }
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
		time=t_taken, kernel=kernel, X_obs=X)
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
