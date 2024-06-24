context("sim_cond_normal")

test_that("`sim_cond_normal` works as expected.",{
  n <- 10
  m <- 2
  joint_mean <- rnorm(n)
  mu_cond <- rnorm(n-m)
  locs_new <- simulatedData$locs[1:m,]
  locs_obs <- simulatedData$locs[1:n,]
  expect_error(sim_cond_normal(joint_mean, mu_cond, locs_new, locs_obs, 
                               kernel=kernel_exp),
               "locs_new and locs_obs must be matrices")
  expect_error(sim_cond_normal(joint_mean, mu_cond, 
                               as.matrix(simulatedData$locs[1:(n+1),]), 
                               as.matrix(locs_obs), 
                               kernel=kernel_exp),
               "Invalid number of rows for `locs_new`.")
  expect_error(sim_cond_normal(joint_mean, mu_cond, 
                               as.matrix(locs_new), 
                               as.matrix(locs_obs), 
                               kernel=kernel_exp, sigma=-1, ell=0),
               "sigma and ell need to be positive.")
  expect_error(sim_cond_normal(joint_mean, mu_cond, 
                               as.matrix(locs_new), 
                               as.matrix(locs_obs), 
                               kernel=kernel_matern, sigma=2, kappa=0),
               "sigma and kappa need to be positive")
  expect_error(kernel_exp(sigma=1, ell=1), 
               "x is not provided. Must provide X1 and X2.")
  expect_error(kernel_matern(sigma=1, kappa=1), 
               "x is not provided. Must provide X1 and X2.")
  expect_error(kernel_exp(sigma=1, ell=1, X1=locs_new, X2=locs_obs), 
               "X1 and X2 must be matrices.")
  expect_error(kernel_matern(sigma=1, kappa=1, X1=locs_new, X2=locs_obs), 
               "X1 and X2 must be matrices.")
})