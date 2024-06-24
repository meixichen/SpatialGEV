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
    # Model a
    fit_a <- spatialGEV_fit(
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
    n_rl <- ifelse(rls[1]==0, 0, length(rls))
    test_model_fit(fit_a, n_loc, n_random=1, n_rl=n_rl)
    # Model ab
    fit_ab <- spatialGEV_fit(
      data = y,
      locs = locs,
      random = "ab",
      init_param = list(
        beta_a = beta_a,
        beta_b = beta_b,
        a = rep(0, n_loc),
        log_b = rep(0, n_loc),
        s = 0,
        log_sigma_a = 0,
        log_kappa_a = 0,
        log_sigma_b = 0,
        log_kappa_b = 0
      ),
      reparam_s = "positive",
      kernel = kernel,
      beta_prior = list(
        beta_a=c(0,100),
        beta_b=c(0,10)
      ),
      X_a = matrix(1, nrow=n_loc, ncol=1),
      X_b = matrix(1, nrow=n_loc, ncol=1),
      matern_pc_prior = list(
        matern_a=matern_pc_prior(1e5,0.95,5,0.1),
        matern_b=matern_pc_prior(1e5,0.95,3,0.1)
      ),
      return_levels = rls,
      silent = TRUE
    )
    test_model_fit(fit_ab, n_loc, s_true=logs, n_random=2, n_rl=n_rl)
  }
})

test_that("Test that `spatialGEV_fit` runs without error for Expo kernel.", {
  rls <- 0.9
  # Model a
  fit_a <- spatialGEV_fit(
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
  test_model_fit(fit_a, n_loc, n_random=1, n_rl=length(rls))
  # Model ab
  fit_ab <- spatialGEV_fit(
    data = y,
    locs = locs,
    random = "ab",
    init_param = list(
      beta_a = beta_a,
      beta_b = beta_b,
      a = rep(0, n_loc),
      log_b = rep(0, n_loc),
      s = logs[1],
      log_sigma_a = 1,
      log_ell_a = 5,
      log_sigma_b = 1,
      log_ell_b = 5
    ),
    reparam_s = "positive",
    kernel = "exp",
    beta_prior = list(
      beta_a=c(0,100),
      beta_b=c(0,10)
    ),
    X_a = matrix(1, nrow=n_loc, ncol=1),
    X_b = matrix(1, nrow=n_loc, ncol=1),
    return_levels = rls,
    silent = TRUE
  )
  test_model_fit(fit_ab, n_loc, s_true=logs, n_random=2, n_rl=length(rls))
})

test_that("`spatialGEV_sample` and `spatialGEV_predict` work for model abs",{
  locs <- simulatedData2$locs
  a <- simulatedData2$a
  logb <- simulatedData2$logb
  logs <- simulatedData2$logs
  y <- simulatedData2$y
  n_loc <- length(y)
  # Set a seed here because we would like to only test the functionality
  # of sample and predict, and not all simulated data has enough complexity
  # to need to be fit using model abs. 
  set.seed(345)
  n_test <- 250
  n_train <- n_loc-n_test
  test_ind <- sample(1:n_loc, n_test)
  train_ind <- setdiff(1:n_loc, test_ind)
  locs_test <- locs[test_ind,]
  y_test <- y[test_ind]
  locs_train <- locs[-test_ind,]
  y_train <- y[-test_ind]
  fit_s <- spatialGEV_fit(data = y_train, locs = locs_train, random = "abs",
                          init_param = list(
                            a = rep(60, n_train), 
                            log_b = rep(3,n_train), 
                            s = rep(-2,n_train),
                            beta_a = 60, beta_b = 2, beta_s = -2,
                            log_sigma_a = 0, log_kappa_a = 0,
                            log_sigma_b = -3, log_kappa_b = -1,
                            log_sigma_s = -1, log_kappa_s = -1),
                          reparam_s = "positive", kernel="spde", silent = T)
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
  method <- "maxsmooth"
  for (random in c("a", "ab")){
    for (kernel in  c("exp", "matern")){
      expect_error(
        fit_a <- spatialGEV_fit(
          data = y,
          locs = locs,
          random = random,
          method = method,
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
          silent = TRUE
        ), 
        "For `method = 'maxsmooth'`, only `random = 'abs'` and `kernel = 'spde'` are currently implemented.")
    }
  }
})
