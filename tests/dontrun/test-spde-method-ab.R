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
fit_matern <- spatialGEV_fit(y = y, locs = locs, random = "ab",
                            init_param = list(beta_a=3, beta_b=1,
					      a = rep(3, n_loc), log_b = rep(0, n_loc), s = logs,
                                              log_sigma_a = 1, log_kappa_a = -2,
					      log_sigma_b = 1, log_kappa_b = -2),
                            kernel = "matern", reparam_s = "positive", silent = F)
print(fit_matern)
ab_matern <- summary(fit_matern$report, "random")
a_matern <- ab_matern[1:n_loc,1]
logb_matern <- ab_matern[(n_loc+1):(2*n_loc),1]

cat("Sampling from the Matern model...\n")
sam_matern <- spatialGEV_sample(fit_matern, 500)
summary_matern <- summary(sam_matern)
#----- End: Without SPDE -----------

#----- With SPDE -------------
fit_spde <- spatialGEV_fit(y=y, locs=locs, random="ab",
                           init_param = list(beta_a=3, beta_b=1,
					     a=rep(3, n_loc),log_b=rep(0, n_loc), s=-1,
					     log_sigma_a=1, log_kappa_a=-2, 
					     log_sigma_b=1, log_kappa_b=-2),
                           kernel = "spde", reparam_s = "positive", silent = F)
ab_spde <- summary(fit_spde$report, "random")
n_s <- nrow(ab_spde)/2 # number of triangles in the mesh
meshidxloc <- fit_spde$meshidxloc # indices of the recorded locations in the mesh
a_spde <- ab_spde[meshidxloc, 1]
logb_spde <- ab_spde[(n_s+1):(2*n_s), ][meshidxloc, 1]

cat("Sampling from the SPDE model...\n")
sam_spde <- spatialGEV_sample(fit_spde, 500)
summary_spde <- summary(sam_spde)
#----- End: With SPDE ----------

#----- Plot the results ----------
par(mfrow=c(2,2))
plot(a, a_matern, xlab="True a", ylab="Estimated a (Matern)")
abline(0, 1, lty="dashed", col="blue")
plot(a, a_spde, xlab="True a", ylab="Estimated a (SPDE)")
abline(0, 1, lty="dashed", col="blue")

plot(logb, logb_matern, xlab="True log(b)", ylab="Estimated log(b) (Matern)") 
abline(0, 1, lty="dashed", col="blue")
plot(logb, logb_spde, xlab="True log(b)", ylab="Estimated log(b) (SPDE)") 
abline(0, 1, lty="dashed", col="blue")
