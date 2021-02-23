#' Fit a hierarchical spatial Gev model.
#'
#' @param y Vector of `n` response values. 
#' @param X `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param random Either "a" or "ab". This indicates which GEV parameters are considered as random effects.
#' @param init.param A list of initial parameters of. See details. 
#' @param reparam.s A flag indicating whether the shape parameter is "zero", "unconstrained", constrained to be "negative", or constrained to be "positive".
#' @param sp.thres Thresholding value to create sparse covariance matrix. Any distance value greater than or equal to `sp.thres` will be set to 0. Default is 0, which means not using sparse matrix.
#' @param adfun.only Only output the ADfun constructed using TMB?
#' @param ignore.random Ignore random effect?
#' @param silent Do not show tracing information?
#' @return If `adfun.only=TRUE`, an list given by `MakeADFun()` is output.
#' If `adfun.only=FALSE`, this function outpus a list containing the following:
#' - An adfun object
#' - A fit object given by calling `nlminb()` on the adfun
#' - An object of class `sdreport` from TMB which contains the point estimates, standard error, and precision matrix for the fixed and random effects
#' @details 
#' This function adopts Laplace approximation using TMB model to integrate out the random effects.
#' Note that `b, s, sigma, ell` are estimated on the log scale.
#' 
#' The log-transformed random effects are assumed to follow Gaussian processes with mean 0 and covariance matrix defined by the exponential covariance function:
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
#' The order of parameters in `init.param` must be: a, log_b, log_s, log_sigma_a, log_ell_a, log_sigma_b, log_ell_b.
#' If reparam.s = "negative" or "postive", the initial value of `s` should be that of log(|s|).
#' @export
spatialGEV_fit <- function(y, X, random, init.param, reparam.s, sp.thres=0, adfun.only=FALSE, ignore.random=FALSE, silent=FALSE){
  
  if (length(y) != nrow(X)){
    stop("The length of y must be the same as the number of rows of X.")
  }
  if (!(random %in% c("a", "ab"))){
    stop("Argument random must be either 'a' or 'ab'.")
  } 
  if (!(reparam.s %in% c("zero", "unconstrained", "positive", "negative"))){
    stop("Argument reparam.s must be one of 'zero', 'unconstrained', 'positive', or 'negative'.")
  } 
  
  mod <- paste("model",random, sep="_")
  n_loc <- length(y)
  dd <- as.matrix(stats::dist(X))
  data <- list(n = n_loc, y = y, dd = dd, sp_thres = sp.thres, reparam_s = reparam.s)
  
  if (random == "ab" & !ignore.random){
    random <- c("a", "log_b")
  }
  else if (ignore.random){
    random <- NULL
  }
  compile(paste0(mod, ".cpp")) ##
  dyn.load(dynlib(mod)) ##
  adfun <- TMB::MakeADFun(data = data,
                          parameters = init.param,
                          random = random,
                          DLL = mod, 
                          silent = silent)
  if (adfun.only){
    adfun
  }
  else{
    fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
    report <- TMB::sdreport(adfun, getJointPrecision = TRUE)
    list(adfun=adfun, fit=fit, report=report) 
  }
}
