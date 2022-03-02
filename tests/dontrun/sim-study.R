devtools::load_all()
library(SpatialExtremes)
set.seed(234)

################ Simulate data #####################
x <- y <- seq(0, 10, length=20)
coor <- cbind(x, y)
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
data <- Map(rgev, n=sample(10:30, n_loc, replace=TRUE),
         loc=as.vector(a), scale=as.vector(b), shape=as.vector(s))
data_mat <- matrix(sapply(data, mean), ncol=sqrt(n_loc))
filled.contour(x, y, data_mat, 
	       color.palette = terrain.colors, main="y")

################ Model fitting using exp #######################
X <- expand.grid(x,y)
init_param <- list(a = rep(50, n_loc), log_b = rep(3, n_loc), s = -2,
                   log_sigma_a = 1, log_ell_a = 1, 
                   log_sigma_b = 1, log_ell_b = 1)
fit_e <- spatialGEV_fit(data, X, random = "ab", 
		      init_param = init_param,
		      reparam_s = "positive", kernel="exp")

################ Model fitting using Matern ##############
init_param <- list(a = rep(50, n_loc), log_b = rep(3, n_loc), s = -2,
                   log_sigma_a = 1, log_kappa_a = -1,
                   log_sigma_b = 1, log_kappa_b = -1)
fit_m <- spatialGEV_fit(data, X, random = "ab",
                      init_param = init_param,
                      reparam_s = "positive", kernel="matern")

################ Model fitting using SPDE #################
init_param <- list(a = rep(50, n_loc), log_b = rep(3, n_loc), s = -2,
                   log_sigma_a = 1, log_kappa_a = -1,
                   log_sigma_b = 1, log_kappa_b = -1)
fit_s <- spatialGEV_fit(data, X, random = "ab",
                      init_param = init_param,
                      reparam_s = "positive", kernel="spde",
                      max.edge = c(1,2))
meshidxloc <- fit_s$meshidxloc
a_s <- fit_s$rep$par.random[meshidxloc]
logb_s <- fit_s$rep$par.random[length(fit_s$rep$par.random)/2+meshidxloc]
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


############### Which locations are not estimated well using SPDE? ############
a_bad <- which((a_s - a) > 2.5)
b_bad <- which((exp(logb_s) - b) < -1.5)

par(mar=c(5.1, 4.1, 4.1, 9.1))
plot(X[,1],X[,2], xlab="Longitude", ylab="Latitude",
     main="Which locations are not estimated well by SPDE")
points(X[a_bad,1], X[a_bad,2], 
       col="red", pch=24, cex=2)
points(X[b_bad,1], X[b_bad,2],
       col="blue", pch=25, cex=2)
legend("topright", inset=c(-0.3,0), 
       legend=c("a","b"), 
       pch=c(24,25), col=c("red","blue"), title="Bad estimates")
