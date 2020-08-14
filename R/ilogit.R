#' The inverse logit function
#'
#' @param x A numeric value.
#' @return The inverse logit of `x`.
#' @export
ilogit <- function(x){
  exp(x)/(1+exp(x))
}