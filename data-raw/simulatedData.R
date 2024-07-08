set.seed(123)
#require(evd)
#require(mvtnorm)
#require(mclust)
# 1. Simulate location coordinates
lon <- seq(0, 10, length.out = 20)
lat <- seq(0, 10, length.out = 20)
locs <- expand.grid(x = lon, y = lat)
n_loc <- nrow(locs)

# 2. Simulate a from the log density surface of bivariate normal
mu_a <- c(4,4)
Sig_a <- diag(2,2)
a <- mvtnorm::dmvnorm(x=locs, mean=mu_a, sigma=Sig_a, log=TRUE)
a <- (a + 30)/5
a_mat <- matrix(a, ncol=sqrt(n_loc))

# 3. Simulate log(b) from the log density surface of bivariate normal mixture of 2 components
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
logb <- mclust::dens("VII", data = locs, logarithm = TRUE, parameters = parameters) 
logb <- (logb - max(logb))/15# rescale log(b)
logb_mat <- matrix(logb, ncol=sqrt(n_loc))

# 4. Simulate s
logs <- -2

# 5. Simulate data
y <- Map(rgev, n=sample(10:30, n_loc, replace=TRUE), 
	 loc=a, scale=exp(logb), shape=exp(logs))

simulatedData <- list(locs=locs, a = a, logb = logb, logs=logs, y = y)
save(simulatedData, file="simulatedData.RData")


#-------- Simulate data with covariates ----------
set.seed(123)
cov_mat <- cbind(rep(1, n_loc), rnorm(n_loc, 0, 0.5), rnorm(n_loc, 1, 0.5))
beta_a <- c(mean(a), rnorm(2, 1, 0.6))
a_withcov <- cov_mat %*% beta_a + a

beta_b <- c(mean(logb), rnorm(1, 1, 0.2))
logb_withcov <- cov_mat[,1:2] %*% beta_b + logb

y_withcov <- Map(rgev, n=sample(10:30, n_loc, replace=TRUE), 
                 loc=a_withcov, scale=exp(logb_withcov), shape=exp(logs))
simulatedData_withcov <- list(locs = locs, a = a_withcov, logb = logb_withcov, 
                              logs = logs, y = y_withcov, beta_a = beta_a,
                              beta_b = beta_b, covariates = cov_mat)
save(simulatedData_withcov, file="simulatedData_withcov.RData")
