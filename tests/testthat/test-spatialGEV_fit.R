context("spatialGEV_fit for model ab")

n_loc <- 20
a <- simulatedData$a[1:n_loc]
logb <- simulatedData$logb[1:n_loc]
logs <- simulatedData$logs[1:n_loc]
y <- simulatedData$y[1:n_loc]
locs <- simulatedData$locs[1:n_loc,]
beta_a <- mean(a)
beta_b <- mean(logb)

test_that("Test that `spatialGEV_fit` runs without error for Matern/SPDE.", {
  rls <- c(0.5, 0.9)
  for (kernel in c("matern", "spde")){
    fit <- spatialGEV_fit(
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
    test_model_fit(fit, n_loc, n_random=2, n_rl=length(rls))
  }
})

test_that("Test that `spatialGEV_fit` runs without error for Expo kernel.", {
  rls <- c(0.5, 0.9)
  fit <- spatialGEV_fit(
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
  test_model_fit(fit, n_loc, n_random=2, n_rl=length(rls))
})