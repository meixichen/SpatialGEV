test_model_fit <- function(fit, n_loc, n_random=2, n_rl=2){
  expect_true(fit$pdHess_avail)
  expect_output(print(fit), "Model fitting")
  
  # Check the summary method for spatialGEV_fit
  fit_summ <- summary(fit)
  expect_named(fit_summ, c("fixed", "random", "return_levels"))
  expect_equal(dim(fit_summ$return_levels), c(n_loc, n_rl*2)) 
  expect_equal(dim(fit_summ$random), c(n_loc*n_random, 2))
  expect_lt(abs((fit_summ$fixed[1,1]-logs[1])/logs[1]), 0.1) # rel. error of s 
  
  #---------- Test the sampling method --------------
  n_draw <- 10
  loc_ind <- sample(n_loc, n_draw, replace = FALSE)
  sam <- spatialGEV_sample(model=fit, n_draw=n_draw,
                           observation=TRUE, loc_ind=loc_ind)
  expect_output(print(sam), "The samples")
  # Check the summary method for spatialGEV_sample
  sam_summ <- summary(sam)
  expect_named(sam_summ, c("param_summary", "y_summary"))
  expect_equal(dim(sam_summ$y_summary), c(n_draw, 6))
  
  #---------- Test the predict method ---------------
  locs_test <- simulatedData$locs[n_loc+1,]
  pred <- spatialGEV_predict(model = fit, locs_new = as.matrix(locs_test),
                             n_draw = n_draw)
  expect_output(print(pred), "posterior predictive")
  # Check the summary method for spatialGEV_predict
  pred_summ <- summary(pred)
  expect_true(all(!is.na(pred_summ)))
  expect_equal(colnames(pred_summ), 
               c("2.5%","25%","50%","75%","97.5%","mean"))
  
  # Predict with pre-run samples
  sam2 <- spatialGEV_sample(model=fit, n_draw=n_draw, loc_ind=1:n_loc)
  pred2 <- spatialGEV_predict(model = fit, locs_new = as.matrix(locs_test),
                              n_draw = n_draw, 
                              parameter_draws = sam2$parameter_draws)
  expect_output(print(pred2), "posterior predictive")
  pred_summ2 <- summary(pred2)
  expect_true(all(!is.na(pred_summ2)))
  expect_equal(colnames(pred_summ2), 
               c("2.5%","25%","50%","75%","97.5%","mean"))
}