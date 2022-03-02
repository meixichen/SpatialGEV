test_that("spatialGEV_fit/sample/predict works fine", {
  # Extract data
  locs <- simulatedData$locs
  a <- simulatedData$a
  logb <- simulatedData$logb
  logs <- simulatedData$logs
  y <- simulatedData$y
  n_loc <- length(y)

  # Test if fitting works
  cat("Fitting...\n")
  fit <- spatialGEV_fit(y = y, locs = locs, random = "a",
			init_param = list(beta_a = mean(a),
					  a = rep(2, n_loc), log_b = -1, s = logs,
					  log_sigma_a = 1, log_ell_a = 5),
                        kernel="exp", X_a=matrix(1,nrow=n_loc,ncol=1),
			reparam_s = "positive", silent = F)
  print(fit)
  expect_equal(fit$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling...\n")
  samps <- spatialGEV_sample(model=fit, n_draw=50, observation=T)
  print(samps)
  samps_summary <- summary(samps)

  # Test if prediction using Matern works
  cat("Predicting...\n")
  n_test <- 20                                # number of test locations
  test_ind <- sample(1:400, n_test)           # indices of the test locations
  locs_test <- locs[test_ind,]                # coordinates of the test locations
  y_test <- y[test_ind]                       # observations at the test locations
  locs_train <- locs[-test_ind,]              # coordinates of the training (observed) locations
  y_train <- y[-test_ind]                     # observations at the training locations
    
  train_fit <- spatialGEV_fit(y = y_train, locs = locs_train, random = "a",
			      init_param = list(beta_a = 3,
						a = rep(2, length(y_train)), log_b = -1, s = logs,
						log_sigma_a = 1, log_kappa_a = -2),
			      kernel="matern", reparam_s = "positive", silent = F)
  pred <- spatialGEV_predict(model = train_fit, locs_new = as.matrix(locs_test),
                             n_draw = 50)
  print(pred)
  pred_summary <- summary(pred)
  success <- 0
  expect_equal(success, 0) 
})
