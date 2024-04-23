context("model_abs_matern")

test_that("`model_abs_matern` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
  n_tests <- 30 # number of test simulations
  for (ii in 1:n_tests){
    # simulate parameters and data
    sim_res <- test_sim(random = c("a", "b", "s"), kernel = "matern", reparam_s = "positive")
    nll_tmb <- calc_tmb_nll(sim_res)
    expect_equal(sim_res$nll_r, nll_tmb)
    # Unconstrained s
    sim_res <- test_sim(random = c("a", "b", "s"), kernel = "matern", reparam_s = "unconstrained")
    nll_tmb <- calc_tmb_nll(sim_res)
    expect_equal(sim_res$nll_r, nll_tmb)
    # Negative s
    sim_res <- test_sim(random = c("a", "b", "s"), kernel = "matern", reparam_s = "negative")
    nll_tmb <- calc_tmb_nll(sim_res)
    expect_equal(sim_res$nll_r, nll_tmb)
  }
})
