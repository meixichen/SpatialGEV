\donttest{
library(SpatialGEV)
n_loc <- 20
a <- simulatedData$a[1:n_loc]
logb <- simulatedData$logb[1:n_loc]
logs <- simulatedData$logs[1:n_loc]
y <- simulatedData$y[1:n_loc]
locs <- simulatedData$locs[1:n_loc,]
beta_a <- mean(a); beta_b <- mean(logb)
fit <- spatialGEV_fit(y, locs = locs, random = "ab",
                      init_param = list(beta_a = beta_a,
                                        beta_b = beta_b,
                                        a = rep(0, n_loc),
                                        log_b = rep(0, n_loc),
                                        s = 0,
                                        log_sigma_a = 0,
                                        log_kappa_a = 0,
                                        log_sigma_b = 0,
                                        log_kappa_b = 0),
                      reparam_s = "positive",
                      kernel = "matern",
                      silent = TRUE)
sam <- spatialGEV_sample(model = fit, n_draw = 100,
                         observation = TRUE, loc_ind=1:10)
print(sam)
summary(sam)
}
