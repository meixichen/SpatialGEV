test_that("spatialGEV_fit/sample/predict works fine", {
  # Extract data
  all_locs <- simulatedData$locs
  smaller_area <- which(all_locs[,1]<5 & all_locs[,2]<5)    
  locs <- all_locs[smaller_area,]
  a <- simulatedData$a[smaller_area]
  logb <- simulatedData$logb[smaller_area]
  logs <- simulatedData$logs
  y <- simulatedData$y[smaller_area]
  n_loc <- length(y)

  n_test <- 10                                # number of test locations
  n_train <- n_loc-n_test
  test_ind <- sample(1:n_loc, n_test)           # indices of the test locations
  locs_test <- locs[test_ind,]                # coordinates of the test locations
  y_test <- y[test_ind]                       # observations at the test locations
  locs_train <- locs[-test_ind,]              # coordinates of the training (observed) locations
  y_train <- y[-test_ind]                     # observations at the training locations
  
  #---------- Test exponential kernel ------------
  
  ##### Model a
  cat("Fitting model a using exp kernel...\n")
  fit_e <- spatialGEV_fit(y = y_train, locs = locs_train, random = "a",
			init_param = list(beta_a = 3,
					  a = rep(2, n_train), log_b = -1, s = logs,
					  log_sigma_a = 1, log_ell_a = 5),
                        kernel="exp",
			reparam_s = "positive", silent = F)
  expect_equal(fit_e$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model a...\n")
  samps_e <- spatialGEV_sample(model=fit_e, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model a...\n")
  pred_e <- spatialGEV_predict(model = fit_e, locs_new = as.matrix(locs_test), n_draw = 5)
  
  ##### Model ab
  cat("Fitting model ab using exp kernel...\n")
  fit_e <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
			init_param = list(beta_a = 3, beta_b = 0,
					  a = rep(2, n_train), log_b = rep(-1,n_train), s = logs,
					  log_sigma_a = 1, log_ell_a = 5,
					  log_sigma_b = 1, log_ell_b = 5),
                        kernel="exp",
			reparam_s = "positive", silent = F)
  expect_equal(fit_e$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model ab...\n")
  samps_e <- spatialGEV_sample(model=fit_e, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model ab...\n")
  pred_e <- spatialGEV_predict(model = fit_e, locs_new = as.matrix(locs_test), n_draw = 5)

  #----------- Test SPDE kernel -----------------
  ###### Model a
  cat("Fitting model a using spde kernel...\n")
  fit_s <- spatialGEV_fit(y = y_train, locs = locs_train, random = "a",
			init_param = list(beta_a = 3,
					  a = rep(2, n_train), log_b = -1, s = logs,
					  log_sigma_a = 1, log_kappa_a = -2),
                        kernel="spde",
			reparam_s = "positive", silent = F)
  expect_equal(fit_s$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model a...\n")
  samps_s <- spatialGEV_sample(model=fit_s, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model a...\n")
  pred_s <- spatialGEV_predict(model = fit_s, locs_new = as.matrix(locs_test), n_draw = 5)
  
  ###### Model ab
  cat("Fitting model ab using SPDE...\n")
  fit_s <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
                        init_param = list(beta_a = 3, beta_b = 0,
                                          a = rep(2, n_train), log_b = rep(-1, n_train), s = logs,
                                          log_sigma_a = 1, log_kappa_a = -2,
					  log_sigma_b = 1, log_kappa_b = -2),
                        kernel="spde",
                        reparam_s = "positive", silent = F)
  print(fit_s)
  expect_equal(fit_s$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model ab...\n")
  samps_s <- spatialGEV_sample(model=fit_s, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model ab...\n")
  pred_s <- spatialGEV_predict(model = fit_s, locs_new = as.matrix(locs_test), n_draw = 5)

  #--------- Test Matern  -----------------------
  ###### Model a
  cat("Fitting model a to the training set using Matern...\n")  
  fit_m <- spatialGEV_fit(y = y_train, locs = locs_train, random = "a",
			      init_param = list(beta_a = 3,
						a = rep(2, n_train), log_b = -1, s = logs,
						log_sigma_a = 1, log_kappa_a = -2),
			      kernel="matern", reparam_s = "positive", silent = F)
  # Test if sampling works
  cat("Sampling from model a...\n")
  samps_m <- spatialGEV_sample(model=fit_m, n_draw=5, observation=T)
  
  cat("Sampling from posterior predictive distribution of model a...\n")
  pred_m <- spatialGEV_predict(model = fit_m, locs_new = as.matrix(locs_test), n_draw = 5)
  
  ###### Model ab
  cat("Fitting model ab using SPDE...\n")
  fit_m <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
                        init_param = list(beta_a = 3, beta_b = 0,
                                          a = rep(2, n_train), log_b = rep(-1, n_train), s = logs,
                                          log_sigma_a = 1, log_kappa_a = -2,
					  log_sigma_b = 1, log_kappa_b = -2),
                        kernel="matern",
                        reparam_s = "positive", silent = F)
  print(fit_m)
  expect_equal(fit_m$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model ab...\n")
  samps_m <- spatialGEV_sample(model=fit_m, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model ab...\n")
  pred_m <- spatialGEV_predict(model = fit_m, locs_new = as.matrix(locs_test), n_draw = 5)
  
  success <- 0
  expect_equal(success, 0) 
})
