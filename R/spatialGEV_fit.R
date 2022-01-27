#' Fit a hierarchical spatial Gev model.
#'
#' @param y List of `n` locations each with `n_obs[i]` independent GEV realizations. 
#' @param X `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param random Either "a" or "ab". This indicates which GEV parameters are considered as random effects.
#' @param init_param A list of initial parameters of. See details. 
#' @param reparam_s A flag indicating whether the shape parameter is "zero", "unconstrained", constrained to be "negative", or constrained to be "positive". See details.
#' @param s_prior Optional. A length 2 vector where the first element is the mean of the normal prior on s or log(s) and the second is the standard deviation.
#' @param kernel Kernel function for spatial random effects covariance matrix. Can be "exp" (exponential kernel), "matern" (Matern kernel), or "spde" (Matern kernel with SPDE approximation).
#' @param sp_thres Optional. Thresholding value to create sparse covariance matrix. Any distance value greater than or equal to `sp_thres` will be set to 0. Default is -1, which means not using sparse matrix. 
#' @param adfun_only Only output the ADfun constructed using TMB?
#' @param ignore_random Ignore random effect?
#' @param silent Do not show tracing information?
#' @param ... Arguments to pass to `INLA::inla.mesh.2d()`.
#' @return If `adfun_only=TRUE`, an list given by `MakeADFun()` is output.
#' If `adfun_only=FALSE`, this function outpus a list containing the following:
#' - An adfun object
#' - A fit object given by calling `nlminb()` on the adfun
#' - An object of class `sdreport` from TMB which contains the point estimates, standard error, and precision matrix for the fixed and random effects
#' @details 
#' This function adopts Laplace approximation using TMB model to integrate out the random effects.
#' 
#' The random effects are assumed to follow Gaussian processes with mean 0 and covariance matrix defined by the chosen kernel function. E.g., using the exponential kernel function:
#' ```
#' cov(i,j) = sigma*exp(-|x_i - x_j|^2/ell)
#' ```
#' Therefore, when specifying the initial parameters, below are some parameters that need to specify.
#' Must specify: `a`, `log_b`, `s`, `log_sigma_a`, `log_ell_a`.
#' Optional: `log_sigma_b`, `log_ell_b`.
#' If the scale parameter b is considered a random effect, its corresponding GP hyperparameters `log_sigma_b` and `log_ell_b` need to be specified.
#' For example, if assume both a and log_b~GP, need to specify 
#' ```
#' init.pram=list(a=rep(1,n),log_b=rep(0,n),s=1,log_sigma_a=0,log_ell_a=0, log_sigma_b=0,log_ell_b=0).
#' ```
#' The order of parameters in `init_param` must be: a, log_b, log_s, log_sigma_a, log_ell_a, log_sigma_b, log_ell_b.
#' If the Matern or SPDE kernel is used, two hyperparameters `sigma` and `kappa` are present for each spatial random effect. I.e., need to specify `sigma_a/b` and `kappa_a/b` in `init_param`.
#' If reparam_s = "negative" or "postive", the initial value of `s` should be that of log(|s|).
#' @export
spatialGEV_fit <- function(y, X, random, init_param, reparam_s, s_prior, kernel="exp", sp_thres=-1, adfun_only=FALSE, ignore_random=FALSE, silent=FALSE){
  
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
    if (kernel == "matern") data$nu <- 1
  }
  else if (kernel == "spde"){
    mesh <- INLA::inla.mesh.2d(X, max.edge=2)
    spde <- (INLA::inla.spde2.matern(mesh)$param.inla)[c("M0", "M1", "M2")]
    n_s <- nrow(spde$M0) # number of mesh triangles created by INLA
    meshidxloc <- as.integer(mesh$idx$loc - 1)
    data <- list(model = mod, y = unlist(y), n_obs = n_obs, meshidxloc = meshidxloc, 
		 reparam_s = reparam_s, spde = spde, nu = 1)
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

  if (missing(s_prior)){
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
    adfun
  }
  else{
    start_t <- Sys.time()
    fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
    report <- TMB::sdreport(adfun, getJointPrecision = TRUE)
    t_taken <- as.numeric(difftime(Sys.time(), start_t, unit="secs"))
    out <- list(adfun=adfun, fit=fit, report=report, time=t_taken, kernel=kernel)
    class(out) <- "spatialGEVfit"
    out 
  }
}
