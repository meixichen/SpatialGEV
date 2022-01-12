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
                            reparam_s = "positive", silent = T)
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
cat("Fitting the model with spde approx...\n")
adfun_with <- MakeADFun(data = list(model = "model_a_spde", 
				       y = y, 
				       meshidxloc = meshidxloc,
				       reparam_s = "positive",
				       s_mean = 9999,
				       s_sd = 9999,
				       spde = spde),
                           parameters = list(a = rep(0, n_s), log_b = -1, s = logs, log_kappa = -2),
                           random = "a",
                           DLL = "SpatialGEV_TMBExports",
                           silent = T)
fit_with <- nlminb(adfun_with$par, adfun_with$fn, adfun_with$gr)
rep_with <- sdreport(adfun_with)
time_with <- difftime(Sys.time(), start)
cat("Time taken to fit with spde approx is", time_with, "\n")
a_with <- summary(rep_with, "random")
a_with <- a_with[meshidxloc, ]
#----- End: With SPDE ----------

#----- Plot the results ----------
par(mfrow=c(1,2))
plot(a, a_with[,1])
plot(a, a_without[,1])
