#' Print method for spatialGEVfit
#'
#' @param x Model object of class `spatialGEVfit` returned by `spatialGEV_fit`.
#' @param ... More arguments for `print`.
#' @return Information about the fitted model containing number of fixed/random effects,
#' fitting time, convergence information, etc.
#' @export

print.spatialGEVfit <- function(x, ...){
  # Time info
  time <- x$time
  ifelse(time <= 600,
	 print_t <- paste(time, "seconds"),
	 print_t <- paste(time/60, "minutes"))
  cat("Model fitting took", print_t, "\n")

  # Convergence info
  if(x$fit$converge == 0){
    mes <- "The model has reached relative convergence \n"
    cat(mes)
  } else {
    cat("The model has not converged and the convergence message output by nlminb is: \n",
        x$fit$message)}

  # Kernel info
  cat("The model uses a", x$kernel, "kernel \n")

  # Parameter info
  cat("Number of fixed effects in the model is", length(x$fit$par), "\n")
  cat("Number of random effects in the model is", length(x$report$par.random), "\n")

  # Hessian info
  ifelse(x$pdHess_avail,
	 mes <- "Hessian matrix is positive definite. Use spatialGEV_sample to obtain posterior samples \n",
         mes <- "Hessian matrix is NOT available or NOT positive definite. spatialGEV_sample and spatialGEV_predict cannot be used \n")
  cat(mes)
}

#' Summary method for spatialGEVfit
#'
#' @param object Object of class `spatialGEVfit` returned by `spatialGEV_fit`.
#' @param ... Additional arguments for `summary`. Not used.
#' @return Point estimates and standard errors of fixed effects, random effects,
#' and the return levels (if specified in `spatialGEV_fit()`) returned by TMB.
#' @export

summary.spatialGEVfit <- function(object, ...){
  fixed_summary <- summary(object$report, "fixed")
  random_summary <- summary(object$report, "random")
  if (object$kernel == "spde"){
    random_len <- nrow(random_summary)/length(object$random)
    loc_ind <- object$meshidxloc
    random_output_ind <- c(loc_ind, loc_ind+random_len, loc_ind+random_len*2)
    random_output_ind <- random_output_ind[random_output_ind <= nrow(random_summary)]
    random_summary <- random_summary[random_output_ind,]
  } else{
    loc_ind <- seq_along(object$locs_obs[,1])
  }
  out <- list(fixed = fixed_summary, random = random_summary)
  if (length(object$return_levels) > 0){
    quantile_summary <- cbind(do.call(cbind, object$return_levels),
                              do.call(cbind, object$return_levels_sd))
    rl_names <- names(object$return_levels)
    colnames(quantile_summary) <- as.vector(
      vapply(c("Estimate", "Std.Error"), function(name) paste0(name, rl_names),
	     letters[seq_along(rl_names)]))
    if (object$kernel == "spde") quantile_summary <- quantile_summary[loc_ind,]
    out$return_levels <- quantile_summary
  }
  out
}

#' Print method for spatialGEVsam
#'
#' @param x Object of class `spatialGEVsam` returned by `spatialGEV_sample`.
#' @param ... Additional arguments for `print`.
#' @return Information about the object including dimension and direction to use `summary` on the object.
#' @export

print.spatialGEVsam <- function(x,...){
  # Dimension info
  dim_param <- dim(x[["parameter_draws"]])
  cat("The samples contain", dim_param[1], "draws of", dim_param[2], "parameters \n")
  if ("y_draws" %in% names(x)){
    dim_y <- dim(x[["y_draws"]])
    cat("The samples contain", dim_y[1], "draws of response at", dim_y[2], "locations \n")
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
#' @param x Object of class `spatialGEVpred` returned by `spatialGEV_predict`.
#' @param ... Additional arguments for `print`.
#' @return Information about the prediction.
#' @export

print.spatialGEVpred <- function(x, ...){
  # Dimension info
  dim_pred <- dim(x$pred_y_draws)
  cat(dim_pred[1], "posterior predictive samples have been draw for", dim_pred[2], "test locations\n")
  cat("The number of training locations is", dim(x$locs_obs)[1], "\n")
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



