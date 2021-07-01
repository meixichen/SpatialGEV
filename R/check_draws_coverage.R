#' Posterior draw coverage check
#'
#' @param y A vector of length `n` containing the true observations. We would like to check if the true observation at each location is covered in the posterior draw.
#' @param draws An `n x k` matrix with `n` being the number of draws from the posterior predictive distribution and `k` being the number of locations
#' @param expected_covr A vector of length `m` containing the expected lengths of credible intervals to look at 
#' @param ci.method Either "ETI" (default) or "HDI" for calculating the credible intervals. The former is the highest density interval and the latter is the equal tail interval.
#' @return A vector of `m`. The i-th entry is the proportion of locations at which the true observation is contained 
#' in the credible interval given by the posterior draws, with the length of credible interval given by the i-th entry 
#' in `expected_covr`.
#' @export
check_draws_coverage <- function(y, draws, expected_covr=seq(0.1, 0.99, by = 0.01), ci.method="ETI"){
  is_covered <- lapply(1:length(expected_covr), 
                       function(j){
                         theor_covr <- expected_covr[j]
                         sapply(1:ncol(draws),
                                function(i){
                                  if (ci.method=="ETI"){
                                    wind_ci <- quantile(draws[,i], c((1-theor_covr)/2, 1-(1-theor_covr)/2))
                                    ifelse({y[i] >= wind_ci[1] & y[i] <= wind_ci[2]}, 1, 0)
                                  }
                                  else if (ci.method=="HDI"){
                                    wind_ci <- bayestestR::ci(draws[,i], ci = theor_covr, method = "HDI")
                                    ifelse({y[i] >= wind_ci$CI_low & y[i] <= wind_ci$CI_high}, 1, 0)
                                  }
                                  else{
                                    stop("ci.method must be one of 'ETI' or 'HDI'.")
                                  }
                                })
                       }) 
  sapply(is_covered, mean) 
}
