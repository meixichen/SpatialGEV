#' Calculate the negative loglikelihood of the hierarchical GEV model.
#'
#' @param a Vector of `n` location paramter
#' @param logit.tau Transformed vector of `n` location paramter by `logit.tau = logit((a-exp(log.s)/exp(log.b))/y_min)`, where `y_min` is a value between 0 and expected minimum of y.
#' @param log.b Log-transformed vector of `n` scale parameters with the constraint of `b>0`. Considered as a random effect.
#' @param log.s Log-transformed vector of `n` shape parameters with the constraint of `s>0`. Considered as a random effect.
#' @param log.siga The squared variance parameter (scalar) on the diagonal of the covariance matrix of `a ~ GP` with the constraint of `sigma_a > 0`.
#' @param log.ella The parameter (scalar) that controls the impact of the distance bewteen two locations on the covariance with the constraint of `ell_a > 0`.
#' @param log.sigb The squared variance parameter (scalar) on the diagonal of the covariance matrix of `logb ~ GP` with the constraint of `sigma_b > 0`.
#' @param log.ellb The parameter (scalar) that controls the impact of the distance bewteen two locations on the covariance with the constraint of `ell_b > 0`.
#' @param log.sigs The squared variance parameter (scalar) on the diagonal of the covariance matrix of `logs ~ GP` with the constraint of `sigma_s > 0`.
#' @param log.ells The parameter (scalar) that controls the impact of the distance bewteen two locations on the covariance with the constraint of `ell_s > 0`.
#' @param y_min The scalar used to transform location parameter into `logit.tau`.
#' @param y Vector of `n` responses of mws values.
#' @param dd A `n x n` distance matrix.
#' @return Scalar value of the negative loglikelihood.
#' @export
hierGev_nll <- function(a, logit.tau, log.b, log.s, log.siga, log.ella, log.sigb, log.ellb, log.sigs, log.ells, y_min, y, dd) {
  mvn_nll_a <- 0
  mvn_nll_b <- 0
  mvn_nll_s <- 0
  if (missing(a)){
    a <- ilogit(logit.tau)*y_min + exp(log.s)/exp(log.b) # recover a from logit.tau
  }
  if (!missing(log.siga)){
    cov.a <- loc_cov(dd = dd, sigma = exp(log.siga), ell = exp(log.ella))
    mvn_nll_a <- -dmvnorm(logit.tau, sigma = cov.a, log = TRUE)
  }
  if (!missing(log.sigb)){
    cov.b <- loc_cov(dd = dd, sigma = exp(log.sigb), ell = exp(log.ellb))
    mvn_nll_b <- -dmvnorm(log.b, sigma = cov.b, log = TRUE)
  }
  if (!missing(log.sigs)){
    cov.s <- loc_cov(dd = dd, sigma = exp(log.sigs), ell = exp(log.ells))
    mvn_nll_s <- -dmvnorm(log.s, sigma = cov.s, log = TRUE)
  }
  gev_nll <- -sum(mapply(dgev, x = y, loc = a, scale = exp(log.b), shape = exp(log.s), MoreArgs = list(log = TRUE)))
  gev_nll + mvn_nll_a + mvn_nll_b + mvn_nll_s
}
