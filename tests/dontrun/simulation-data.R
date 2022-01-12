require(evd) # GEV density and simulation
require(mvtnorm)
require(mclust)
set.seed(123)
# 1. Simulate location coordinates
lon <- seq(0, 10, length.out = 20)
lat <- seq(0, 10, length.out = 20)
locs <- expand.grid(x = lon, y = lat)
n_loc <- nrow(locs)

# 2. Simulate log(b) from the log density surface of bivariate normal mixture
G <- 2 # number of components
pro <- c(0.4,0.6) # mixing proportion
mu1 <- c(1, 0) # mean of the first component
mu2 <- c(8, 7) # mean of the second component
mean <- cbind(mu1, mu2)
diag_entry <- c(0.5, 1)
Sig1 <- diag(diag_entry[1], 2) # covariance matrix of the first component
Sig2 <- diag(diag_entry[2], 2) # covariance matrix of the second component
variance <- list(d=2, G=G, sigma=array(Sig1, Sig2), sigmasq=diag_entry, scale=diag_entry)
parameters <- list(pro=pro, mean=mean, variance=variance)
#simdat <- mclust::simVII(parameters=parameters, n=200, seed=123)
# Simulate logb 
logb <- mclust::dens("VII", data = locs, logarithm = TRUE, parameters = parameters) 
logb <- (logb - max(logb))/15# rescale log(b)
logb_mat <- matrix(logb, ncol=sqrt(n_loc))
# persp(x=lon, y=lat, z=logb_mat, theta=45, phi=35, r=10, expand=1, axes=T,
#       ticktype="detailed", xlab="lon", ylab="lat", zlab="log(b)", main="Spatial spread of log(b)")
# filled.contour(x=lon, y=lat, z=logb_mat, key.title = title(main="log(b)"), 
#                main = "Contour plot of log(b)", xlab="lon", ylab="lat")

# 3. Simulate a from the log density surface of bivariate normal
mu_a <- c(4,4)
Sig_a <- diag(2,2)
a <- mvtnorm::dmvnorm(x=locs, mean=mu_a, sigma=Sig_a, log=TRUE)
a <- (a + 30)/5
a_mat <- matrix(a, ncol=sqrt(n_loc))
# persp(x=lon, y=lat, z=a_mat, theta=45, phi=35, r=10, expand=1, axes=T,
#       ticktype="detailed", xlab="lon", ylab="lat", zlab="a", main="Spatial spread of a")
# filled.contour(x=lon, y=lat, z=a_mat, key.title = title(main="a"), 
#                main = "Contour plot of a", xlab="lon", ylab="lat")

# 4. Simulate s
logs <- -2
# 5. Simulate wind
y <- unlist(Map(rgev, n=1, loc=a, scale=exp(logb), shape=exp(logs)))
y_mat <- matrix(y, ncol=sqrt(n_loc))
# persp(x=lon, y=lat, z=y_mat, theta=45, phi=35, r=10, expand=1, axes=T,
#       ticktype="detailed", xlab="lon", ylab="lat", zlab="y", main="Spatial spread of simulated y")
# filled.contour(x=lon, y=lat, z=y_mat, key.title = title(main="y"), 
#                main = "Contour plot of wind", xlab="lon", ylab="lat")
# par(mfrow=c(1,1))
# hist(y)
