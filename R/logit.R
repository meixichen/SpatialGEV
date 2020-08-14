#' The logit function
#'
#' @param x A numeric value.
#' @return The logit-transformed `x`.
#' @export
logit <- function(x){
  log(x/(1-x))
}