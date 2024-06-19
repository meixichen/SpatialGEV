#' @export
#' @importFrom posterior as_draws_df
as_draws_df.spatialGEVsam <- function(x, ...){
  foo <- function(x){
    matches <- str_match(x, '([abs])([0-9]{1,})')
    new <- str_c(matches[,2], "[", matches[,3], "]")
    ifelse(is.na(new), x, new)
  }
  colnames(x$parameter_draws) <- foo(colnames(x$parameter_draws))
  as_draws_df(x$parameter_draws)
}

#' @export
#' @importFrom broom tidy
tidy.spatialGEVfit <- function(x, effects = c("fixed", "random")){
  effects <- match.arg(effects)
  model_summary <- summary(x$report, select = effects)
  ret <- as_tibble(model_summary, rownames = "term")
  names(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")
  if(identical(effects, "random")){
    ret <- ret |>
      mutate(
        term = sprintf("%s[%s]", term, rownames(x$locs_obs))
      )
  }
  ret
}
