\donttest{
library(SpatialGEV)
n_loc <- 20
a <- simulatedData$a[1:n_loc]
logb <- simulatedData$logb[1:n_loc]
logs <- simulatedData$logs[1:n_loc]
y <- simulatedData$y[1:n_loc]
locs <- simulatedData$locs[1:n_loc,]
# No covariates are included, only intercept is included.
fit <- spatialGEV_fit(
  data = y,
  locs = locs,
  random = "ab",
  init_param = list(
    a = rep(0, n_loc),
    log_b = rep(0, n_loc),
    s = 0,
    beta_a = 0,
    beta_b = 0,
    log_sigma_a = 0,
    log_kappa_a = 0,
    log_sigma_b = 0,
    log_kappa_b = 0
  ),
  reparam_s = "positive",
  kernel = "matern",
  X_a = matrix(1, nrow=n_loc, ncol=1),
  X_b = matrix(1, nrow=n_loc, ncol=1),
  silent = TRUE
)
print(fit)

# To use a different optimizer other than the default `nlminb()`, create
# an object ready to be passed to optimizer functions using `adfun_only=TRUE`
obj <- spatialGEV_fit(
  data = y,
  locs = locs, random = "ab",
  init_param = list(
    a = rep(0, n_loc),
    log_b = rep(0, n_loc),
    s = 0,
    beta_a = 0,
    beta_b = 0,
    log_sigma_a = 0,
    log_kappa_a = 0,
    log_sigma_b = 0,
    log_kappa_b = 0
  ),
  reparam_s = "positive",
  kernel = "matern",
  X_a = matrix(1, nrow=n_loc, ncol=1),
  X_b = matrix(1, nrow=n_loc, ncol=1),
  adfun_only = TRUE
)
fit <- optim(obj$par, obj$fn, obj$gr)
}

# Using the SPDE kernel (SPDE approximation to the Matern kernel)
# Make sure the INLA package is installed before using `kernel="spde"`
\dontrun{
library(INLA)
n_loc <- 20
y <- simulatedData2$y[1:n_loc]
locs <- simulatedData2$locs[1:n_loc,]
fit_spde <- spatialGEV_fit(
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
  kernel = "spde",
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
plot(fit_spde$mesh) # Plot the mesh
points(locs[,1], locs[,2], col="red", pch=16) # Plot the locations
}
