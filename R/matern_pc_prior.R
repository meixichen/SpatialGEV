#' Helper funcion to specify a Penalized Complexity (PC) prior on the Matern hyperparameters
#' 
#' @param rho_0 Hyperparameter for PC prior on the range parameter. Must be positive. 
#' See details.
#' @param p_rho Hyperparameter for PC prior on the range parameter. Must be between 0
#' and 1. See details.
#' @param sig_0 Hyperparameter for PC prior on the range parameter. Must be positive. 
#' See details.
#' @param p_sig Hyperparameter for PC prior on the range parameter. Must be between 0
#' and 1. See details.
#' @return A list to provide to the `matern_pc_prior` argument of `spatialGEV_fit`.
#' @details
#' The joint prior on `rho` and `sig` achieves
#' ```
#' P(rho < rho_0) = p_rho,
#' ```
#' and
#' ```
#' P(sig > sig_0) = p_sig,
#' ```
#' where `rho = sqrt(8*nu)/kappa` and `sig = sqrt(sigma)`.
#' @examples
#' \dontrun{
#' fit_spde <- spatialGEV_fit(y = y, locs = locs, random = "abs",
#'                            init_param = list(a = rep(0, n_loc),
#'                                              log_b = rep(0, n_loc), 
#'                                              s = rep(-2, n_loc),
#'                                              beta_a = 0,
#'                                              beta_b = 0,
#'                                              beta_s = -2,
#'                                              log_sigma_a = 0, 
#'                                              log_kappa_a = 0,
#'                                              log_sigma_b = 0, 
#'                                              log_kappa_b = 0,
#'                                              log_sigma_s = 0, 
#'                                              log_kappa_s = 0,
#'                                              ),
#'                            reparam_s = "positive",
#'                            kernel = "spde",
#'                            beta_prior = list(beta_a=c(0,100), beta_b=c(0,10),
#'                                              beta_s=c(0,10)),
#'                            matern_pc_prior = list(
#'                                                  matern_a=matern_pc_prior(1e5,0.95,5,0.1),
#'                                                  matern_b=matern_pc_prior(1e5,0.95,3,0.1),
#'                                                  matern_s=matern_pc_prior(1e2,0.95,1,0.1),
#'                                                  )
#'                            adfun_only = TRUE) 
#' }
#' @export
matern_pc_prior <- function(rho_0, p_rho, sig_0, p_sig){
  if (p_rho > 1 | p_rho < 0) stop("p_rho must be between 0 and 1.")
  if (p_sig > 1 | p_sig < 0) stop("p_sig must be between 0 and 1.")
  if (rho_0 < 0 | sig_0 < 0) stop("rho_0 and sig_0 must be positive.")
  out <- list(range_prior=c(rho_0, p_rho),
              sigma_prior=c(sig_0, p_sig)) 
  class(out) <- "PC_prior"
  return(out)
}
