#require(SpatialGEV)
require(INLA)
require(TMB)
source("simulation-data.R")
#------ Without SPDE ---------
start <- Sys.time()
cat("Fitting the model without spde approx...\n")
fit_without <- spatialGEV_fit(y = y, X = locs, random = "a",
                            init_param = list(a = rep(0, n_loc), log_b = -1, s = logs,
                                              log_sigma_a = 1, log_ell_a = 5),
                            kernel="exp", reparam_s = "positive", silent = T)
time_without <- difftime(Sys.time(), start)
cat("Time taken to fit without spde approx is", time_without, "\n")
a_without <- summary(fit_without$report, "random")
#----- End: Without SPDE -----------

#----- With SPDE -------------
mesh <- inla.mesh.2d(as.matrix(locs), max.edge=2)
spde <- (inla.spde2.matern(mesh)$param.inla)[c("M0", "M1", "M2")]
n_s <- nrow(spde$M0)
meshidxloc <- mesh$idx$loc - 1
start <- Sys.time()
fit_with <- spatialGEV_fit(y=y, X=locs, random="a", 
			   init_param = list(a=rep(0,400),log_b=-1,s=-1,log_kappa=-2), 
			   kernel = "spde", reparam_s = "positive", silent = T)
time_with <- difftime(Sys.time(), start)
cat("Time taken to fit with spde approx is", time_with, "\n")
a_with <- summary(fit_with$report, "random")
a_with <- a_with[meshidxloc, ]
#----- End: With SPDE ----------

#----- Plot the results ----------
par(mfrow=c(1,2))
plot(a, a_with[,1])
plot(a, a_without[,1])
