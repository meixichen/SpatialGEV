#' @export
spatial_gev <- function(x, ...){
  UseMethod("spatial_gev")
}

#' @export
spatial_gev.formula <- function(
    formula, data, coordinates = ~ x + y,
    method = c("laplace", "maxsmooth"),
    kernel = c("spde", "matern", "exp"),
    reparam_s = c("zero", "unconstrained", "negative", "positive"),
    prior = list(),
    control = list(),
    init = NULL
){
  method <- match.arg(method)
  kernel <- match.arg(kernel)
  
  formula <- Formula::Formula(formula)
  formula <- expand_dot(formula)
  flen <- length(formula)
  stopifnot(identical(flen[1], 1L), flen[2] <= 3 & flen[2])
  
  lhs <- rlang::f_lhs(formula)
  data <- tidyr::chop(data, lhs) # make data nloc x nvar with a list response
  
  Xmats <- regression_matrix(formula, data)
  
  ret <- list(
    data = get(lhs, data),
    locs = model.frame(coordinates, data),
    random  = switch(flen[2], "a", "ab", "abc"),
    method = method,
    reparam_s = reparam_s,
    kernel = kernel,
    formula = formula,
    control = set_control(control),
    prior   = set_prior(prior)
  )
  if(is.null(init)){
    ret$init_param <- init_params(
      x = ret$data,
      random = Filter(Negate(is.null), lapply(Xmats, ncol)),
      kernel = ret$kernel,
      use_fgev = TRUE
    )
  }
  ret <- c(ret, setNames(Xmats, paste0("X_", names(Xmats))))
  mycall <- rlang::call2("spatialGEV_fit", !!!ret[c(1:6,9:12)], !!!ret$control, !!!ret$prior, .ns = "SpatialGEV")
  fit <- eval(mycall)
  attr(fit, "formula") <- ret$formula
  attr(fit, "call") <- mycall
  fit
}

set_prior <- function(prior){
  ret <- list(
    nu              = 1,
    s_prior         = NULL,
    beta_prior      = NULL,
    matern_pc_prior = NULL
  )
  nms <- names(ret)
  nms_set <- names(prior)
  ret[nms_set] <- prior
  if (length(noNms <- nms_set[!nms_set %in% nms])){
    stop("Unknown control parameters passed: ", paste(noNms, collapse = ", "))
  }
  ret
}

set_control <- function(control){
  ret <- list(
    return_levels         = 0,
    get_return_levels_cov = TRUE,
    sp_thres              = -1,
    adfun_only            = FALSE,
    ignore_random         = FALSE,
    get_hessian           = TRUE,
    silent                = TRUE,
    mesh_extra_init = list(
      a     = 0,
      log_b = -1,
      s     = 0.001
    )
  )
  nms <- names(ret)
  nms_set <- names(control)
  ret[nms_set] <- control
  if (length(noNms <- nms_set[!nms_set %in% nms])){
    stop("Unknown control parameters passed: ", paste(noNms, collapse = ", "))
  }
  ret
}

#' @export
init_params <- function(x, random, kernel, use_fgev = TRUE){
  # For each parameter need:
  # list("a", "log_b", "s") where length is 1 when not random and nloc when random
  # For each random parameter above, need
  hyper <- if(identical(kernel, "exp")){
    c("log_sigma", "log_ell")
  } else {
    c("log_sigma", "log_kappa")
  }
  random_nms <- names(random)
  betas <- paste0("beta_", random_nms)
  hyper <- paste0(hyper, "_", rep(random_nms, each = 2))
  
  n_loc <- length(x)
  x <- unlist(x)
  if(use_fgev){
    start <- evd::fgev(x)
    start <- unname(mvtnorm::rmvnorm(1, start$estimate, vcov(start))[1, , drop = TRUE])
  } else {
    m <- mean(x)
    s <- sqrt(6 * var(x)) / pi
    start <- c(m - 0.57722 * s, s, 1e-8)
  }
  
  params <- c(
    list(
      a     = rep(start[1], n_loc),
      log_b = log(start[2]),
      s     = start[3]
    ),
    as.list(setNames(rep(0, length(betas)), betas)),
    as.list(setNames(rep(0, length(hyper)), hyper))
  )
  params$beta_a <- rep(params$beta_a, random$a)
  if(is.element("b", random_nms)){
    params$log_b <- rep(params$log_b, n_loc)
    params$beta_b <- rep(params$beta_b, random$b)
  }
  if(is.element("s", random_nms)){
    params$s <- rep(params$s, n_loc)
    params$beta_s <- rep(params$beta_s, random$s)
  }
  params
}

# Formula parsing helpsers -----------------------------------------------------

expand_dot <- function(x){
  ret <- attr(terms(x, dot = "previous"), "Formula_without_dot")
  ret <- if(is.null(ret)) x else ret
  ret
}

#' Extract design matrix from all RHS random effect terms
#'
#' Extract the design matrix from the data.frame for all RHS random effect terms.
#' For the spatial GEV model here, any RHS terms should consist solely of a single
#' RE term, optionally with regression terms.
regression_matrix <- function(x, data){
  f <- \(i){
    ret <- formula(x, lhs = 0, rhs = i)
    ret <- get_design_formula(ret)[[1]] #Only allow one RE term per rhs
    model.matrix(ret, data)
  }
  ret <- setNames(rep(list(NULL), 3), c('a', 'b', 's'))
  n <- length(x)[2]
  ret[1:n] <- lapply(1:n, f)
  ret
}

#' Extract the regression matrix component of a lme style random effect term
#' i.e., the dots in ( ... | group). Returns a list of formulas with as many values as there
#' are RE terms in the formula.
get_design_formula <- function(x){
  bar <- lme4::findbars(x)
  design <- lapply(bar, rlang::call_args)
  lapply(design, \(x) eval(call("~", x[[1]])))
}

