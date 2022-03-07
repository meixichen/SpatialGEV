test_that("spatialGEV_fit/sample/predict works fine for model a and ab", {
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
			init_param = list(a = rep(2, n_train), log_b = -1, s = logs,
					  beta_a =3, log_sigma_a = 1, log_ell_a = 5),
                        kernel="exp",
			reparam_s = "positive", silent = T)
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
			init_param = list(a = rep(2, n_train), log_b = rep(-1,n_train), s = logs,
					  beta_a = 3, beta_b = 0,
					  log_sigma_a = 1, log_ell_a = 5,
					  log_sigma_b = 1, log_ell_b = 5),
                        kernel="exp",
			reparam_s = "positive", silent = T)
  expect_equal(fit_e$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model ab...\n")
  samps_e <- spatialGEV_sample(model=fit_e, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model ab...\n")
  pred_e <- spatialGEV_predict(model = fit_e, locs_new = as.matrix(locs_test), n_draw = 5)
  
  ##### Model abs
  cat("Fitting model abs using exp kernel...\n")
  fit_e <- spatialGEV_fit(y = y_train, locs = locs_train, random = "abs",
			init_param = list(
					  a = rep(2, n_train), 
					  log_b = rep(-1,n_train), 
					  s = rep(-3,n_train),
					  beta_a = 3, beta_b = 0, beta_s = -2,
					  log_sigma_a = 1, log_ell_a = 5,
					  log_sigma_b = 1, log_ell_b = 5,
					  log_sigma_s = -1, log_ell_s = 5),
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
			init_param = list(
					  a = rep(2, n_train), log_b = -1, s = logs,
					  beta_a = 3,
					  log_sigma_a = 1, log_kappa_a = -2),
                        kernel="spde",
			reparam_s = "positive", silent = T)
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
                        init_param = list(
                                          a = rep(2, n_train), log_b = rep(-1, n_train), s = logs,
					  beta_a = 3, beta_b = 0,
                                          log_sigma_a = 1, log_kappa_a = -2,
					  log_sigma_b = 1, log_kappa_b = -2),
                        kernel="spde",
                        reparam_s = "positive", silent = T)
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
			      init_param = list(
						a = rep(2, n_train), log_b = -1, s = logs,
						beta_a = 3,
						log_sigma_a = 1, log_kappa_a = -2),
			      kernel="matern", reparam_s = "positive", silent = T)
  # Test if sampling works
  cat("Sampling from model a...\n")
  samps_m <- spatialGEV_sample(model=fit_m, n_draw=5, observation=T)
  
  cat("Sampling from posterior predictive distribution of model a...\n")
  pred_m <- spatialGEV_predict(model = fit_m, locs_new = as.matrix(locs_test), n_draw = 5)
  
  ###### Model ab
  cat("Fitting model ab using SPDE...\n")
  fit_m <- spatialGEV_fit(y = y_train, locs = locs_train, random = "ab",
                        init_param = list(
                                          a = rep(2, n_train), log_b = rep(-1, n_train), s = logs,
					  beta_a = 3, beta_b = 0,
                                          log_sigma_a = 1, log_kappa_a = -2,
					  log_sigma_b = 0, log_kappa_b = -2),
                        kernel="matern",
                        reparam_s = "positive", silent = T)
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


test_that("spatialGEV_fit/sample/predict works fine for model abs", {
  # Extract data
  locs <- simulatedData2$locs
  a <- simulatedData2$a
  logb <- simulatedData2$logb
  logs <- simulatedData2$logs
  y <- simulatedData2$y
  n_loc <- length(y)

  n_test <- 100                                # number of test locations
  n_train <- n_loc-n_test
  test_ind <- sample(1:n_loc, n_test)           # indices of the test locations
  train_ind <- setdiff(1:n_loc, test_ind)
  locs_test <- locs[test_ind,]                # coordinates of the test locations
  y_test <- y[test_ind]                       # observations at the test locations
  locs_train <- locs[-test_ind,]              # coordinates of the training (observed) locations
  y_train <- y[-test_ind]                     # observations at the training locations
  
  #---------- Test exponential kernel ------------
  
  cat("Fitting model abs using exp kernel...\n")
  fit_e <- spatialGEV_fit(y = y_train, locs = locs_train, random = "abs",
			  init_param = list(
					    a = rep(60, n_train), 
					    log_b = rep(2,n_train), 
					    s = rep(-3,n_train),
					    beta_a = 60, beta_b = 2, beta_s = -2,
					    log_sigma_a = 1.5, log_ell_a = 5,
					    log_sigma_b = 1.5, log_ell_b = 5,
					    log_sigma_s = -1, log_ell_s = 5),
			  reparam_s = "positive", kernel="exp", silent = T)
  expect_equal(fit_e$fit$convergence, 0)

  # Test if sampling works
  cat("Sampling from model abs...\n")
  samps_e <- spatialGEV_sample(model=fit_e, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model abs...\n")
  pred_e <- spatialGEV_predict(model = fit_e, locs_new = as.matrix(locs_test), 
			       n_draw = 5, parameter_draws = samps_e$parameter_draws)

  #----------- Test SPDE kernel -----------------

  cat("Fitting model abs using SPDE kernel...\n")
  fit_s <- spatialGEV_fit(y = y_train, locs = locs_train, random = "abs",
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
  cat("Sampling from model abs...\n")
  samps_s <- spatialGEV_sample(model=fit_s, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model abs...\n")
  pred_s <- spatialGEV_predict(model = fit_s, locs_new = as.matrix(locs_test),
			       n_draw = 5, parameter_draws = samps_s$parameter_draws)

  #--------- Test Matern  -----------------------

  cat("Fitting model abs using Matern kernel...\n")
  fit_m <- spatialGEV_fit(y = y_train, locs = locs_train, random = "abs",
			  init_param = list(
					    a = rep(60, n_train), 
					    log_b = rep(3,n_train), 
					    s = rep(-2,n_train),
					    beta_a = 60, beta_b = 3, beta_s = -2,
					    log_sigma_a = 0, log_kappa_a = 0,
					    log_sigma_b = -5, log_kappa_b = -1,
					    log_sigma_s = -1, log_kappa_s = 0),
			  reparam_s = "positive", kernel="matern", silent = T)
  expect_equal(fit_m$fit$convergence, 0)
  
  # Test if sampling works
  cat("Sampling from model abs...\n")
  samps_m <- spatialGEV_sample(model=fit_m, n_draw=5, observation=T)

  # Test if prediction works
  cat("sampling from posterior predictive distribution of model abs...\n")
  pred_m <- spatialGEV_predict(model = fit_m, locs_new = as.matrix(locs_test), 
			       n_draw = 5, parameter_draws = samps_m$parameter_draws)
  
  success <- 0
  expect_equal(success, 0) 
})
