#' Print method for spatialGEVfit
#'
#' @param object Model object of class `spatialGEVfit` returned by `spatialGEV_fit`.
#' @param ... More arguments for `print`.
#' @return Information about the fitted model containing number of fixed/random effects, 
#' fitting time, convergence information, etc.
#' @export

print.spatialGEVfit <- function(x, ...){
  # Time info
  time <- object$time
  ifelse(time <= 600,
	 print_t <- paste(time, "seconds"),
	 print_t <- paste(time/60, "minutes"))
  cat("Model fitting took", print_t, "\n")

  # Convergence info
  if(object$fit$converge == 0){
    mes <- "The model has reached relative convergence \n"
    cat(mes)
  } else {
    cat("The model has not converged and the convergence message output by nlminb is: \n",
        object$fit$message)}
  
  # Kernel info
  cat("The model uses a", object$kernel, "kernel \n")

  # Parameter info
  cat("Number of fixed effects in the model is", length(object$fit$par), "\n")
  cat("Number of random effects in the model is", length(object$report$par.random), "\n")

  # Hessian info
  ifelse(object$report$pdHess, 
	 mes <- "Hessian matrix is positive definite. Use spatialGEV_sample to obtain posterior samples \n",
         mes <- "Hessian matrix is NOT positive definite. spatialGEV_sample and spatialGEV_predict cannot be used \n")
  cat(mes)
}

#' Print method for spatialGEVsam
#' 
#' @param object Object of class `spatialGEVsam` returned by `spatialGEV_sample`.
#' @param ... Additional arguments for `print`.
#' @return Information about the object including dimension and direction to use `summary` on the object.
#' @export 

print.spatialGEVsam <- function(x,...){
  # Dimension info
  dim_param <- dim(object[["parameter_draws"]])
  cat("The samples contains", dim_param[1], "draws of", dim_param[2], "parameters \n")
  if ("y_draws" %in% names(object)){
    dim_y <- dim(object[["y_draws"]])
    cat("The samples contains", dim_y[1], "draws of response at", dim_y[2], "locations \n")
  }
  cat("Use summary() to obtain summary statistics of the samples \n")
}

#' Summary method for spatialGEVsam
#'
#' @param object Object of class `spatialGEVsam` returned by `spatialGEV_sample`.
#' @param q A vector of quantile values used to summarize the samples.
#' Default is `c(0.025, 0.25, 0.5, 0.75, 0.975)`.
#' @param ... Additional arguments for `summary`. Not used.
#' @return Summary statistics of the posterior samples.
#' @export

summary.spatialGEVsam <- function(object, q=c(0.025, 0.25, 0.5, 0.75, 0.975), ...){
  # Summary of all parameters
  param_samps <- object[["parameter_draws"]]
  param_summary <- t(apply(param_samps, 2, quantile, 
	             probs=q))
  param_summary <- cbind(param_summary, apply(param_samps, 2, mean))
  colnames(param_summary)[ncol(param_summary)] <- "mean"
  # Summary of data resamples 
  if ("y_draws" %in% names(object)){
    y_samps <- object[["y_draws"]]
    y_summary <- t(apply(y_samps, 2, quantile, 
	           probs=q))
    y_summary <- cbind(y_summary, apply(y_samps, 2, mean))
    colnames(y_summary)[ncol(y_summary)] <- "mean"
    out <- list(param_summary=param_summary, y_summary=y_summary)
  }
  else{
    out <- list(param_summary=param_summary)
  }
  out
}

#' Print method for spatialGEVpred
#'
#' @param object Object of class `spatialGEVpred` returned by `spatialGEV_predict`.
#' @param ... Additional arguments for `print`.
#' @return Information about the prediction.
#' @export

print.spatialGEVpred <- function(x, ...){
  # Dimension info
  dim_pred <- dim(object$pred_y_draws)
  cat(dim_pred[1], "posterior predictive samples have been draw for", dim_pred[2], "test locations\n")
  cat("The number of training locations is", dim(object$locs_obs)[1], "\n")
  cat("Use summary() to obtain summary statistics of the posterior predictive samples \n")
}

#' Summary method for spatialGEVpred
#'
#' @param object Object of class `spatialGEVpred` returned by `spatialGEV_predict`.
#' @param q A vector of quantile values used to summarize the samples.
#' Default is `c(0.025, 0.25, 0.5, 0.75, 0.975)`.
#' @param ... Additional arguments for `summary`.
#' @return Summary statistics of the posterior predictive samples.
#' @export

summary.spatialGEVpred <- function(object, q=c(0.025, 0.25, 0.5, 0.75, 0.975), ...){
  # Summary of all parameters
  y_draws <- object[["pred_y_draws"]]
  out <- t(apply(y_draws, 2, quantile, 
	             probs=q))
  out <- cbind(out, apply(y_draws, 2, mean))
  colnames(out)[ncol(out)] <- "mean"
  out
}



