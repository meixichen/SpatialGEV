#' Helper funcion to specify a Penalized Complexity (PC) prior on the Matern hyperparameters
#'
#' @param rho_0 Hyperparameter for PC prior on the range parameter. Must be positive.
#' See details.
#' @param p_rho Hyperparameter for PC prior on the range parameter. Must be between 0
#' and 1. See details.
#' @param sig_0 Hyperparameter for PC prior on the scale parameter. Must be positive.
#' See details.
#' @param p_sig Hyperparameter for PC prior on the scale parameter. Must be between 0
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
#' where `rho = sqrt(8*nu)/kappa`.
#' @references
#' Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017).
#' Penalising model component complexity: A principled, practical approach to
#' construct priors. Statistical Science.
#' @example examples/matern_pc_prior.R
#' @export
matern_pc_prior <- function(rho_0, p_rho, sig_0, p_sig) {
  if (any(vapply(list(rho_0, p_rho, sig_0, p_sig), length, 1L) != 1L)){
    stop("The inputs must be scalars.")
  }
  if (p_rho > 1 | p_rho < 0) stop("p_rho must be between 0 and 1.")
  if (p_sig > 1 | p_sig < 0) stop("p_sig must be between 0 and 1.")
  if (rho_0 < 0 | sig_0 < 0) stop("rho_0 and sig_0 must be positive.")
  out <- list(range_prior=c(rho_0, p_rho),
              sigma_prior=c(sig_0, p_sig))
  class(out) <- "PC_prior"
  return(out)
}
