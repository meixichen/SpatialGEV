#require(SpatialGEV)
require(INLA)
require(TMB)
source("simulation-data.R")
#------ Without SPDE ---------
start <- Sys.time()
cat("Fitting the model without spde approx...\n")
fit_matern <- spatialGEV_fit(y = y, X = locs, random = "a",
                            init_param = list(a = rep(2, n_loc), log_b = -1, s = logs,
                                              log_sigma_a = 1, log_kappa_a = -2),
                            kernel="matern", reparam_s = "positive", silent = F)
time_matern <- difftime(Sys.time(), start)
cat("Time taken to fit without spde approx is", time_matern, "\n")
a_matern <- summary(fit_matern$report, "random")
#----- End: Without SPDE -----------

#----- With SPDE -------------
start <- Sys.time()
cat("Fitting the model with spde approx...\n")
fit_spde <- spatialGEV_fit(y=y, X=locs, random="a",
                           init_param = list(a=rep(2, n_loc),log_b=-1, s=-1,
                                             log_sigma_a=1, log_kappa_a=-2),
                           kernel = "spde", reparam_s = "positive", silent = F)
time_spde <- difftime(Sys.time(), start)
cat("Time taken to fit with spde approx is", time_spde, "\n")
a_spde <- summary(fit_spde$report, "random")
meshidxloc <- fit_spde$adfun$env$data$meshidxloc # indices of the recorded locations in the mesh
a_spde <- a_spde[meshidxloc, ]
#----- End: With SPDE ----------

#----- SPDE with quadrature --------
start <- Sys.time()
cat("Fitting the model with spde approx and quadrature...\n")
adfun_quad <- spatialGEV_fit(y=y, X=locs, random="a",
                           init_param = list(a=rep(2, n_loc),log_b=-1, s=-1,
                                             log_sigma_a=1, log_kappa_a=-2),
                           kernel = "spde", reparam_s = "positive", adfun_only = T, silent = F)
fit_quad <- aghq::marginal_laplace_tmb(adfun_quad, 3, adfun_quad$par) # use aghq's quadrature method
samples_quad <- aghq::sample_marginal(fit_quad, 1000)
time_quad <- difftime(Sys.time(), start)
cat("Time taken to fit with spde approx and quadrature is", time_quad, "\n")
a_samps <- samples_quad$samps
a_quad <- cbind(apply(a_samps, 1, mean)[meshidxloc], apply(a_samps, 1, sd)[meshidxloc])

#----- End: with quadrature

#----- Plot the results ----------
par(mfrow=c(1,3))
plot(a, a_matern[1:n_loc,1], xlab="True a", ylab="Estimated a (Matern)")
abline(0, 1, lty="dashed", col="blue")

plot(a, a_spde[,1], xlab="True a", ylab="Estimated a (SPDE)")
abline(0, 1, lty="dashed", col="blue")

plot(a, a_quad[,1], xlab="True a", ylab="Estimated a (SPDE w/ quad)")
abline(0, 1, lty="dashed", col="blue")

