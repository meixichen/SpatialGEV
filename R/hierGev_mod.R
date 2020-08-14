#' Construct a TMB model object for the spatial hierarchical model
#'
#' @param y Vector of `n` response values. 
#' @param x `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param theta A list of initial parameters of. See details.
#' @param random A vector of parameter names that are considered as random effect. Has to be one or multiple of "a", "b", "s".
#' @param sp.thres Thresholding value to create sparse covariance matrix. Any distance value greater than or equal to `sp.thres` will be set to 0. Default is 0, which means not using sparse matrix.
#' @return A TMB model object, which is a list with components (fn, gr, etc) suitable for calling an R optimizer
#' @details 
#' This function is very similar to `hierGev_fit` except that it doesn't perform optimization inside the function.
#' Therefore, it can be used when the user are interested in exploring the TMB model instead of getting the fits through a black box.
#' Note that `b, s, sigma, ell` are estimated on the log scale.
#' To avoid likelihood of GEV bing 0, instead of directly using the location parameter `a` in estimation, its transformed version is used
#' ```
#' logit_tau = logit((a-s/b)/y_min)
#' ```
#' where `y_min` is any number between `0` and the minimum value of `y`.
#' 
#' The log-transformed random effects are assumed to follow Gaussian processes with mean 0 and covariance matrix defined by the exponential covariance function:
#' ```
#' cov(i,j) = sigma*exp(-|x_i - x_j|^2/ell)
#' ```
#' Therefore, when specifying the initial parameters, below are some parameters that need to specify.
#' Must specify: `logit_tau`, `log_b`, `log_s`.
#' Optional: `log_sigma_a`, `log_sigma_b`, `log_sigma_s`, `log_ell_a`, `log_ell_b`, `log_ell_s`.
#' If one of the above is considered as following GP, its corresponding `sigma` and `ell` need to be specified.
#' For example, if assume log_b~GP, need to specify `theta=c(logit_tau=1,log_b=rep(0,n),log_s=1,log_sigma_b=0,log_ell_b=0)`.
#' @export
hierGev_mod <- function(y, x, theta, random, sp.thres=0){
  
  candidates1 <- c("a","b","s")
  inx <- which(candidates1 %in% random)
  chosen <- candidates1[inx]
  mod <- paste("model",paste(chosen,collapse=""), sep="_")
  n <- length(y)
  dd <- as.matrix(dist(x))
  data <- list(model=mod, n=n, y=y, dd=dd, sp_thres=sp.thres)
  
  if ("a" %in% random){
    candidates2 <- c("logit_tau", "log_b", "log_s")
    random <- candidates2[inx]
    y_min <- runif(1,0, min(y))
    data[["y_min"]] <- y_min
    model <- TMB::MakeADFun(data=data,
                            parameters=theta,
                            random=random,
                            DLL = "SpatialGEV_TMBExports", 
                            silent = TRUE)
  }
  else {
    candidates2 <- c("a", "log_b", "log_s")
    random <- candidates2[inx]
    model <- TMB::MakeADFun(data=data,
                            parameters=theta,
                            random=random,
                            DLL = "SpatialGEV_TMBExports", 
                            silent = TRUE)
  }
  model
}
