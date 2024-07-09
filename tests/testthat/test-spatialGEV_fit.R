context("spatialGEV_fit")

n_loc <- 20
a <- simulatedData$a[1:n_loc]
logb <- simulatedData$logb[1:n_loc]
logs <- simulatedData$logs[1:n_loc]
y <- simulatedData$y[1:n_loc]
locs <- simulatedData$locs[1:n_loc,]
beta_a <- mean(a)
beta_b <- mean(logb)
rl_list <- list(0, c(0.5, 0.9))
kernel_list <- c("matern", "spde")
test_settings <- expand.grid(rl=rl_list, kernel=kernel_list)[-3,]

test_that("Test that `spatialGEV_fit` runs without error for Matern/SPDE.", {
  for (ii in seq_along(test_settings[,1])){
    kernel <- as.character(test_settings$kernel)[ii]
    rls <- test_settings$rl[[ii]]
    fit_args <- list(
      data = y,
      locs = locs,
      random = "a",
      init_param = list(
        beta_a = beta_a,
        a = rep(0, n_loc),
        log_b = 0,
        s = 0,
        log_sigma_a = 0,
        log_kappa_a = 0
      ),
      reparam_s = "zero",
      kernel = kernel,
      return_levels = rls,
      silent = TRUE
    )
    # Model a
    fit_a <- do.call(spatialGEV_fit, fit_args)
    n_rl <- ifelse(rls[1]==0, 0, length(rls))
    test_model_fit(fit_a, n_loc, n_random=1, n_rl=n_rl)
    # Model ab
    fit_args$random <- "ab"
    fit_args$init_param <- list(beta_a = beta_a, beta_b = beta_b,
                                a = rep(0, n_loc), log_b = rep(0, n_loc),  
                                s = 0,
                                log_sigma_a = 0, log_kappa_a = 0,
                                log_sigma_b = 0, log_kappa_b = 0)
    fit_args$reparam_s <- "positive"
    fit_args$beta_prior <- list(beta_a=c(0,100), beta_b=c(0,10))
    fit_args$X_a <- fit_args$X_b <- matrix(1, nrow=n_loc, ncol=1)
    fit_args$matern_pc_prior <- list(matern_a=matern_pc_prior(1e5,0.95,5,0.1),
                                     matern_b=matern_pc_prior(1e5,0.95,3,0.1))
    fit_ab <- do.call(spatialGEV_fit, fit_args)
    test_model_fit(fit_ab, n_loc, s_true=logs, n_random=2, n_rl=n_rl)
  }
})

test_that("Test that `spatialGEV_fit` runs without error for Expo kernel.", {
  rls <- 0.9
  # Model a
  fit_args <- list(
    data = y,
    locs = locs,
    random = "a",
    init_param = list(
      beta_a = beta_a,
      a = rep(0, n_loc),
      log_b = 0,
      s = 0,
      log_sigma_a = 1,
      log_ell_a = 5
    ),
    reparam_s = "zero",
    kernel = "exp",
    beta_prior = list(
      beta_a=c(0,100)
    ),
    return_levels = rls,
    silent = TRUE
  )
  fit_a <- do.call(spatialGEV_fit, fit_args)
  test_model_fit(fit_a, n_loc, n_random=1, n_rl=length(rls))
  # Model ab
  fit_args$random <- "ab"
  fit_args$init_param <- list(beta_a = beta_a, beta_b = beta_b, 
                              a = rep(0, n_loc), log_b = rep(0, n_loc),
                              s = logs[1],
                              log_sigma_a = 1, log_ell_a = 5,
                              log_sigma_b = 1, log_ell_b = 5)
  fit_args$reparam_s <- "positive"
  fit_args$beta_prior <- list(beta_a=c(0,100), beta_b=c(0,10))
  fit_ab <- do.call(spatialGEV_fit, fit_args)
  test_model_fit(fit_ab, n_loc, s_true=logs, n_random=2, n_rl=length(rls))
})

test_that("`spatialGEV_sample` and `spatialGEV_predict` work for model abs",{
  reparam_s <- "positive"
  locs <- simulatedData3$locs
  a <- simulatedData3$a
  logb <- simulatedData3$logb
  logs <- simulatedData3$logs
  y <- simulatedData3$y
  n_loc <- length(y)
  n_test <- 2
  n_train <- n_loc-n_test
  test_ind <- sample(1:n_loc, n_test)
  train_ind <- setdiff(1:n_loc, test_ind)
  locs_test <- locs[test_ind,]
  y_test <- y[test_ind]
  locs_train <- locs[-test_ind,]
  y_train <- y[-test_ind]
  fit_args <- list(data = y_train, locs = locs_train, random = "abs",
                   init_param = list(
                     a = a[train_ind], 
                     log_b = logb[train_ind], 
                     s = logs[train_ind],
                     beta_a = mean(a), beta_b = mean(logb), beta_s = mean(logs),
                     log_sigma_a = 0, log_ell_a = 0,
                     log_sigma_b = -3, log_ell_b = 0,
                     log_sigma_s = -1, log_ell_s = 0),
                   reparam_s = reparam_s, kernel="exp", silent=TRUE)
  #--------------- Exponential kernel -------------------
  fit_e <- do.call(spatialGEV_fit, fit_args)
  # Test if prediction works
  pred_e <- spatialGEV_predict(model = fit_e, 
                               locs_new = as.matrix(locs_test[1:2,]),
                               n_draw = 2)
  expect_true(all(!is.na(summary(pred_e))))
  
  #------------ SPDE kernel --------------
  fit_args$init_param <- c(fit_args$init_param[1:6],
                           log_sigma_a = 0, log_kappa_a = 0,
                           log_sigma_b = -3, log_kappa_b = -1,
                           log_sigma_s = -1, log_kappa_s = -1)
  fit_args$kernel <- "spde"
  fit_args$max.edge <- 2
  fit_s <- do.call(spatialGEV_fit, fit_args)
  expect_equal(fit_s$fit$convergence, 0)
  # Test if sampling works
  samps_s <- spatialGEV_sample(model=fit_s, n_draw=2)
  expect_named(summary(samps_s), "param_summary")
  # Test if prediction works
  pred_s <- spatialGEV_predict(model = fit_s, 
                               locs_new = as.matrix(locs_test[1:2,]),
                               n_draw = 2, 
                               parameter_draws = samps_s$parameter_draws)
  expect_true(all(!is.na(summary(pred_s))))
})


test_that("`spatialGEV_fit` gives an informative error message for features not 
          implemented.", {
  # Only implemented model abs and kernel spde for max-and-smooth method
  init_param <- list(
    beta_a = beta_a,
    a = rep(0, n_loc),
    log_b = 0,
    s = 0,
    log_sigma_a = 0,
    log_kappa_a = 0
  )
  for (random in c("a", "ab")){
    for (kernel in  c("exp", "matern")){
      expect_error(
        spatialGEV_fit(
          data = y,
          locs = locs,
          random = random,
          method = "maxsmooth",
          init_param = init_param,
          reparam_s = "zero",
          kernel = kernel,
          silent = TRUE
        ), 
        "For `method = 'maxsmooth'`, only `random = 'abs'` and `kernel = 'spde'` are currently implemented.")
    }
  }
  #----- Test for informative errors for maxsmooth -----#
  fit_args <- list(data = list(est=NULL),
                   locs = locs,
                   random = "abs",
                   method = "maxsmooth",
                   init_param = init_param,
                   reparam_s = "zero",
                   kernel = "spde",
                   silent = TRUE)
  expect_error(
    do.call(spatialGEV_fit, fit_args), 
    "For `method == 'maxsmooth'`, `data` must be a list with named elements 'est' and 'var'.")
  
  fit_args$data <- list(est=matrix(rnorm(4), 2, 2), var=matrix(rnorm(4), 2, 2))
  expect_error(
    do.call(spatialGEV_fit, fit_args), 
    "Incorrect dimensions for `data$est`.", fixed=TRUE)
  
  fit_args$data <- list(est=matrix(rnorm(n_loc*3), n_loc, 3), 
                        var=matrix(rnorm(4), 2, 2))
  expect_error(
    do.call(spatialGEV_fit, fit_args), 
    "Incorrect dimensions for `data$var`.", fixed=TRUE)
  
  #----- Test for informative errors for the laplace method -----#
  fit_args$data <- list(rep(NA, 2))
  fit_args$method <- "laplace"
  expect_error(
    do.call(spatialGEV_fit, fit_args), 
    "For `method == 'laplace', `data must be a numeric list.")
  
  fit_args$data <- list(rnorm(1))
  expect_error(
    do.call(spatialGEV_fit, fit_args), 
    "For `method == 'laplace', must have `length(data) == nrow(locs)`.", fixed=TRUE)  
  
  #----- Test for informative errors for arguments -----#
  fit_args$data <- y
  fit_args$reparam_s <- "constrained"
  expect_error(do.call(spatialGEV_fit, fit_args),
               "Argument reparam_s must be one of 'zero', 'unconstrained', 'positive', or 'negative'.")
  
  fit_args$reparam_s <- "zero"
  expect_error(do.call(spatialGEV_fit, fit_args),
               "When s is a random effect, reparam_s cannot be zero.")
  
  fit_args$beta_prior <- c(0, 10)
  fit_args$reparam_s <- "positive"
  expect_error(do.call(spatialGEV_fit, fit_args), "Check beta_prior.")
  
  fit_args$beta_prior <- NULL
  fit_args$matern_pc_prior <- c(0.1, 99)
  expect_error(do.call(spatialGEV_fit, fit_args),
               "Check matern_pc_prior: must be a named list*")
  
  fit_args$matern_pc_prior <- list(matern_a = list(c(1000, 0.95), c(5, 0.1)))
  expect_error(do.call(spatialGEV_fit, fit_args),
               "`matern_pc_prior$matern_a` must be created with `matern_pc_prior()`.",
               fixed=TRUE)
})


test_that("The max-smooth method 
          (a feature for publication purposes and hence not maintained) 
          is working properly.", {
  # Note: the max-smooth feature only works for model abs with kernel spde and
  # positive s parameter. This feature will not be further developed as it was
  # for publication purposes only.
  n_loc <- 50
  locs <- simulatedData3$locs[1:n_loc,]
  a <- simulatedData3$a[1:n_loc]
  logb <- simulatedData3$logb[1:n_loc]
  logs <- simulatedData3$logs[1:n_loc]
  y <- simulatedData3$y[1:n_loc]
  fit_s <- spatialGEV_fit(data = y, method="maxsmooth",
                          locs = locs, random = "abs",
                          init_param = list(
                            a = a, 
                            log_b = logb, 
                            s = logs,
                            beta_a = mean(a), beta_b = mean(logb), beta_s = mean(logs),
                            log_sigma_a = 0, log_kappa_a = 0,
                            log_sigma_b = -3, log_kappa_b = -1,
                            log_sigma_s = -1, log_kappa_s = -1),
                          s_prior=c(0,10),
                          reparam_s = "positive", kernel="spde", silent = T)
  expect_true(fit_s$pdHess_avail)
})
