#' The negative loglikelihood of the hierarchical GEV model given by TMB
#'
#' @param y Vector of `n` responses of mws values.
#' @param dd A `n x n` distance matrix.
#' @param y_min The scalar used to transform location parameter into `logit.tau`.
#' @param theta A list of parameters to be evaluated at in the TMB model.
#' @param random A vector of parameter names that are considered as random effect. Has to be one or multiple of "a", "b", "s".
#' @return Scalar value of the negative loglikelihood.
#' @export
hierGev_tmb_nll <- function(y, dd, y_min, theta, random){
  candidates1 <- c("a","b","s")
  inx <- which(candidates1 %in% random)
  chosen <- candidates1[inx]
  mod <- paste("model",paste(chosen,collapse=""), sep="_")
  
  data <- list(model=mod, n=length(y), y=y, dd=dd, sp_thres=0)
  
  if ("a" %in% random){
    data[["y_min"]] <- y_min
  }
  
  model <- TMB::MakeADFun(data=data,
                     parameters=theta,
                     DLL = "SpatialGEV_TMBExports", 
                     silent = TRUE)
  model$fn(unlist(theta))
}
