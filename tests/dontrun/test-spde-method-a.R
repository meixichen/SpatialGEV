#require(SpatialGEV)
require(INLA)
require(TMB)
locs <- simulatedData$locs
a <- simulatedData$a
logb <- simulatedData$logb
logs <- simulatedData$logs
y <- simulatedData$y
n_loc <- length(y)
#------ Without SPDE ---------
cat("Fitting the model without spde approx...\n")
fit_matern <- spatialGEV_fit(y = y, locs = locs, random = "a",
                            init_param = list(beta_a=3, a = rep(2, n_loc), log_b = -1, s = logs,
                                              log_sigma_a = 1, log_kappa_a = -2),
                            kernel="matern", reparam_s = "positive", silent = F)
print(fit_matern)
a_matern <- summary(fit_matern$report, "random")[1:n_loc,1]
#----- End: Without SPDE -----------

#----- With SPDE -------------
cat("Fitting the model with spde approx...\n")
fit_spde <- spatialGEV_fit(y=y, locs = locs, random="a",
                           init_param = list(beta_a=3, a=rep(2, n_loc),log_b=-1, s=-1,
                                             log_sigma_a=1, log_kappa_a=-2),
                           kernel = "spde", reparam_s = "positive", silent = F)
print(fit_spde)
a_spde <- summary(fit_spde$report, "random")
meshidxloc <- fit_spde$meshidxloc # indices of the recorded locations in the mesh
a_spde <- a_spde[meshidxloc, 1]
#----- End: With SPDE ----------

#----- Plot the results ----------
par(mfrow=c(1,2))
plot(a, a_matern, xlab="True a", ylab="Estimated a (Matern)")
abline(0, 1, lty="dashed", col="blue")

plot(a, a_spde, xlab="True a", ylab="Estimated a (SPDE)")
abline(0, 1, lty="dashed", col="blue")

