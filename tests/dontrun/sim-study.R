devtools::load_all()
library(SpatialExtremes)
set.seed(234)

################ Simulate data #####################
x <- y <- seq(0, 10, length=30)
coor <- cbind(x, y)
locs <- expand.grid(x,y)
n_loc <- length(x)*length(y)
# Simulate a
## nugget is a constant offset, sill is sigma, range is 2*ell 
a <- rgp(1, coor, cov.mod="powexp", 
	  mean=60, nugget=0, sill=5, range=10, smooth=1,
	  grid=TRUE)
filled.contour(x, y, a, color.palette = terrain.colors, main="a")

# Simulate b
b <- rgp(1, coor, cov.mod="powexp", 
          mean=20, nugget=0, sill=3, range=10, smooth=1,
          grid=TRUE)
filled.contour(x, y, b, color.palette = terrain.colors, main="b")


# Simulate s 
s <- rgp(1, coor, cov.mod="powexp",
          mean=0.1, nugget=0, sill=0.001, range=10, smooth=1,
          grid=TRUE)
filled.contour(x, y, s, color.palette = terrain.colors, main="s")

# Simulate data
data <- Map(rgev, n=sample(30:50, n_loc, replace=TRUE),
         loc=as.vector(a), scale=as.vector(b), shape=as.vector(s))
data_mat <- matrix(sapply(data, mean), ncol=sqrt(n_loc))
filled.contour(x, y, data_mat, 
	       color.palette = terrain.colors, main="y")

################ Model fitting using exp #######################
init_param <- list(
		   a = rep(50, n_loc), 
		   log_b = rep(3, n_loc), 
		   s = rep(-2,n_loc),
		   beta_a = 50, beta_b = 20, beta_s = -2, 
                   log_sigma_a = 1, log_ell_a = 1, 
                   log_sigma_b = 1, log_ell_b = 1,
                   log_sigma_s = -1, log_ell_s = 1)
fit_e <- spatialGEV_fit(data, locs, random = "abs", 
		      init_param = init_param,
		      reparam_s = "positive", kernel="exp")
saveRDS(fit_e, "fit_e900.rds")
a_e <- fit_e$rep$par.random[1:n_loc]
logb_e <- fit_e$rep$par.random[(n_loc+1):(2*n_loc)]
logs_e <- fit_e$rep$par.random[(n_loc*2+1):(3*n_loc)]
################ Model fitting using Matern ##############
init_param <- list(
		   a = rep(50, n_loc), 
		   log_b = rep(3, n_loc), 
		   s = rep(-2,n_loc),
		   beta_a = 50, beta_b = 20, beta_s = -2, 
                   log_sigma_a = 1, log_kappa_a = -1,
                   log_sigma_b = 1, log_kappa_b = -1,
                   log_sigma_s = 1, log_kappa_s = -1)
fit_m <- spatialGEV_fit(data, locs, random = "abs",
                      init_param = init_param,
                      reparam_s = "positive", kernel="matern")
saveRDS(fit_m, "fit_m900.rds")
a_m <- fit_m$rep$par.random[1:n_loc]
logb_m <- fit_m$rep$par.random[(n_loc+1):(2*n_loc)]
logs_m <- fit_m$rep$par.random[(n_loc*2+1):(3*n_loc)]

################ Model fitting using SPDE #################
init_param <- list(
		   a = rep(50, n_loc), 
		   log_b = rep(3, n_loc), 
		   s = rep(-2,n_loc),
		   beta_a = 50, beta_b = 20, beta_s = -2, 
                   log_sigma_a = 1, log_kappa_a = -1,
                   log_sigma_b = 1, log_kappa_b = -1,
                   log_sigma_s = -1, log_kappa_s = -1)
fit_s <- spatialGEV_fit(data, locs, random = "abs",
                      init_param = init_param,
                      reparam_s = "positive", kernel="spde")
saveRDS(fit_s, "fit_s900.rds")
meshidxloc <- fit_s$meshidxloc
a_s <- fit_s$rep$par.random[meshidxloc]
logb_s <- fit_s$rep$par.random[length(fit_s$rep$par.random)/3+meshidxloc]
logs_s <- fit_s$rep$par.random[length(fit_s$rep$par.random)/3*2+meshidxloc]
################ Plotting #############################
par(mfrow=c(3,2))
plot(as.vector(a), fit_e$rep$par.random[1:n_loc],
     xlab="True a", ylab="Estimated a", main="True vs estimated a (Exp)")
abline(0, 1, lty="dashed", col="blue")
plot(as.vector(b), exp(fit_e$rep$par.random[(1+n_loc):(2*n_loc)]),
     xlab="True b", ylab="Estimated b", main="True vs estimated b (Exp)")
abline(0, 1, lty="dashed", col="blue")

plot(as.vector(a), fit_m$rep$par.random[1:n_loc],
     xlab="True a", ylab="Estimated a", main="True vs estimated a (Matern)")
abline(0, 1, lty="dashed", col="blue")
plot(as.vector(b), exp(fit_m$rep$par.random[(1+n_loc):(2*n_loc)]),
     xlab="True b", ylab="Estimated b", main="True vs estimated b (Matern)")
abline(0, 1, lty="dashed", col="blue")

plot(as.vector(a), a_s,
     xlab="True a", ylab="Estimated a", main="True vs estimated a (SPDE)")
abline(0, 1, lty="dashed", col="blue")
plot(as.vector(b), exp(logb_s),
     xlab="True b", ylab="Estimated b", main="True vs estimated b (SPDE)")
abline(0, 1, lty="dashed", col="blue")

# Time plot
par(mfrow=c(1,1))
barplot(c(fit_e$time, fit_m$time, fit_s$time), 
        names.arg=c("Exponential", "Matern", "Matern with SPDE approx"), 
	ylab="Time (sec)")


