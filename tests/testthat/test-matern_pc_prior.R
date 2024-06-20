context("matern_pc_prior")

test_that("`matern_pc_prior` gives results of the correct format.",{
  for (ii in seq_len(10)){
    prior_list <- matern_pc_prior(runif(1, 1, 1e5),runif(1, 0, 1),
                                  runif(1, 0, 1e3),runif(1, 0, 1))
    expect_equal(class(prior_list), "PC_prior")
    expect_true(is.list(prior_list))
    expect_true(all(vapply(prior_list, function(x){all(!is.na(x))}, TRUE)))
  }
})

test_that("`matern_pc_prior` throw errors when given out-of-range inputs.",{
  for (ii in seq_len(10)){
    wrong_rho_0 <- runif(1, -10, 0)
    wrong_p_rho <- runif(1, 1, 10)
    wrong_sig_0 <- runif(1, -10, 0)
    wrong_p_sig <- runif(1, -10, 0)
    expect_error(matern_pc_prior(runif(1, 1, 1e5),wrong_p_rho, 
                                 runif(1, 0, 1e3),runif(1, 0, 1)),
                 "p_rho must be between 0 and 1.")
    expect_error(matern_pc_prior(runif(1, 1, 1e5),runif(1, 0, 1), 
                                 runif(1, 0, 1e3),wrong_p_sig),
                 "p_sig must be between 0 and 1.")
    expect_error(matern_pc_prior(wrong_rho_0,runif(1, 0, 1), 
                                 runif(1, 0, 1e3),runif(1, 0, 1)),
                 "rho_0 and sig_0 must be positive.")
    expect_error(matern_pc_prior(runif(1, 1, 1e5),runif(1, 0, 1), 
                                 wrong_sig_0,runif(1, 0, 1)),
                 "rho_0 and sig_0 must be positive.")
    expect_error(matern_pc_prior(sample(1:10, 2),sample(1:10, 2),
                                 sample(1:10, 2),sample(1:10, 2)),
                 "The inputs must be scalars.")
  }
})
