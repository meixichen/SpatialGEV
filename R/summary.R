#' Print method for spatialGEVfit
#'
#' @param object Model object of class `spatialGEVfit` returned by `spatialGEV_fit`.
#' @param ... More arguments for `print`.
#' @return Information about the fitted model containing number of fixed/random effects, 
#' fitting time, convergence information, etc.
#' @export

print.spatialGEVfit <- function(object, ...){
  # Time info
  time <- object$time
  ifelse(time <= 600,
	 print_t <- paste(time, "seconds."),
	 print_t <- paste(time/60, "minutes."))
  cat("Model fitting took", print_t, "\n")

  # Convergence info
  if(object$fit$converge == 0){
    mes <- "The model has reached relative convergence. \n"
    cat(mes)
  } else {
    cat("The model has not converged and the convergence message output by nlminb is: \n",
        object$fit$message)}

  # Parameter info
  cat("Number of fixed effects in the model is", length(object$fit$par), ". \n")
  cat("Number of random effects in the model is", length(object$report$par.random), ". \n")

  # Hessian info
  ifelse(object$report$pdHess, 
	 mes <- "Hessian matrix is positive definite. Use spatialGEV_sample to obtain posterior samples. \n",
         mes <- "Hessian matrix is NOT positive definite. spatialGEV_sample and spatialGEV_predict cannot be used.")
  cat(mes)
}


#' Summary method for spatialGEVsam
#'
#' @param object Object of class `spatialGEVsam` returned by `spatialGEV_sample`.
#' @param ... Additional arguments for `summary`.
#' @return Summary statistics of the posterior samples.
#' @export

summary.spatialGEVsam <- function(object, ...){
  # Summary of all parameters
  param_samps <- object[["parameter_draws"]]
  param_summary <- t(apply(param_samps, 2, quantile, 
	             probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  param_summary <- cbind(param_summary, apply(param_samps, 2, mean))
  colnames(param_summary)[ncol(param_summary)] <- "mean"
  # Summary of data resamples 
  if ("y_draws" %in% names(object)){
    y_samps <- object[["y_draws"]]
    y_summary <- t(apply(y_samps, 2, quantile, 
	           probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
    y_summary <- cbind(y_summary, apply(y_samps, 2, mean))
    colnames(y_summary)[ncol(y_summary)] <- "mean"
    out <- list(param_summary=param_summary, y_summary=y_summary)
  }
  else{
    out <- list(param_summary=param_summary)
  }
  out
}
