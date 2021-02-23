#' Calculate the negative loglikelihood of the spatial GEV model.
#'
#' @param y Vector of `n` independent GEV realizations.
#' @param dd An `n x n` distance matrix.
#' @param a Vector of `n` location paramter
#' @param log_b A numeric value or a vector of `n` log-transformed scale parameters if considered as a random effect.
#' @param s Shape parameter
#' @param log_sigma_a The squared variance parameter (scalar) on the diagonal of the covariance matrix of `a ~ GP` with the constraint of `sigma_a > 0`.
#' @param log_ell_a The parameter (scalar) that controls the impact of the distance bewteen two locations on the covariance with the constraint of `ell_a > 0`.
#' @param log_sigma_b The squared variance parameter (scalar) on the diagonal of the covariance matrix of `logb ~ GP` with the constraint of `sigma_b > 0`.
#' @param log_ell_b The parameter (scalar) that controls the impact of the distance bewteen two locations on the covariance with the constraint of `ell_b > 0`.
#' @return Scalar value of the negative loglikelihood.
#' @details This function is used to test if TMB and R output the same negative loglikelihood.
#' @export
r_nll <- function(y, dd, a, log_b, s, log_sigma_a, log_ell_a, log_sigma_b, log_ell_b) {
  cov_a <- exp(log_sigma_a)*exp(-dd/exp(log_ell_a))
  nll <- -dmvnorm(a, sigma = cov_a, log = TRUE) 
  if (length(log_b)==1){
    nll <- nll - sum(mapply(dgev, x = y, loc = a, 
                            MoreArgs = list(scale = exp(log_b), shape = s, log = TRUE)))
  }
  else { # if b is considered a random effect
    cov_b <- exp(log_sigma_b)*exp(-dd/exp(log_ell_b))
    nll <- nll - dmvnorm(log_b, sigma = cov_b, log = TRUE)
    nll <- nll - sum(mapply(dgev, x = y, loc = a, scale = exp(log_b), 
                 MoreArgs = list(shape = s, log = TRUE)))
  }
  nll
}
