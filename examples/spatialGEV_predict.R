\donttest{
set.seed(123)
library(SpatialGEV)
n_loc <- 20
a <- simulatedData$a[1:n_loc]
logb <- simulatedData$logb[1:n_loc]
logs <- simulatedData$logs[1:n_loc]
y <- simulatedData$y[1:n_loc]
locs <- simulatedData$locs[1:n_loc,]
n_test <- 5
test_ind <- sample(1:n_loc, n_test)

# Obtain coordinate matrices and data lists
locs_test <- locs[test_ind,]
y_test <- y[test_ind]
locs_train <- locs[-test_ind,]
y_train <- y[-test_ind]

# Fit the GEV-GP model to the training set
train_fit <- spatialGEV_fit(
  data = y_train,
  locs = locs_train,
  random = "ab",
  init_param = list(
    beta_a = mean(a),
    beta_b = mean(logb),
    a = rep(0, n_loc-n_test),
    log_b = rep(0, n_loc-n_test),
    s = 0,
    log_sigma_a = 1,
    log_kappa_a = -2,
    log_sigma_b = 1,
    log_kappa_b = -2
  ),
  reparam_s = "positive",
  kernel = "matern",
  silent = TRUE
)

pred <- spatialGEV_predict(
  model = train_fit,
  locs_new = locs_test,
  n_draw = 100
)
summary(pred)
}
