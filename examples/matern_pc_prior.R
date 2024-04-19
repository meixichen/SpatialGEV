\donttest{
n_loc <- 20
y <- simulatedData2$y[1:n_loc]
locs <- simulatedData2$locs[1:n_loc,]
fit <- spatialGEV_fit(
  data = y,
  locs = locs,
  random = "abs",
  init_param = list(
    a = rep(0, n_loc),
    log_b = rep(0, n_loc),
    s = rep(-2, n_loc),
    beta_a = 0,
    beta_b = 0,
    beta_s = -2,
    log_sigma_a = 0,
    log_kappa_a = 0,
    log_sigma_b = 0,
    log_kappa_b = 0,
    log_sigma_s = 0,
    log_kappa_s = 0
  ),
  reparam_s = "positive",
  kernel = "matern",
  beta_prior = list(
    beta_a=c(0,100),
    beta_b=c(0,10),
    beta_s=c(0,10)
  ),
  matern_pc_prior = list(
    matern_a=matern_pc_prior(1e5,0.95,5,0.1),
    matern_b=matern_pc_prior(1e5,0.95,3,0.1),
    matern_s=matern_pc_prior(1e2,0.95,1,0.1)
  )
)
}
