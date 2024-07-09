set.seed(123)
# 1. Simulate location coordinates
lon <- seq(0, 5, length.out = 10)
lat <- seq(0, 5, length.out = 10)
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
logb <- mclust::dens("VII", data = locs, logarithm = TRUE, parameters = parameters)
logb <- (logb +40)/12# rescale log(b)
b <- exp(logb)
logb_mat <- matrix(logb, ncol=sqrt(n_loc))

# 3. Simulate a from the log density surface of bivariate normal
mu_a <- c(3,3)
Sig_a <- diag(3,2)
a <- mvtnorm::dmvnorm(x=locs, mean=mu_a, sigma=Sig_a, log=TRUE)
a <- (a + 80)
a_mat <- matrix(a, ncol=sqrt(n_loc))

# 4. Simulate s
mu_s <- c(4,4)
Sig_s <- diag(2,2)
logs <- mvtnorm::dmvnorm(x=locs, mean=mu_s, sigma=Sig_s, log=TRUE)
logs <- (logs)/5
s <- exp(logs)
logs_mat <- matrix(logs, ncol=sqrt(n_loc))

# 5. Simulate observed data
data <- Map(rgev, n=sample(20:50, n_loc, replace=T),
            loc=a, scale=exp(logb), shape=s)

# 6. 90% quantile of data at each location
z_true <- unlist(Map(qgev, p=0.1, loc=a, scale=b,
                     shape=s, lower.tail=F))

simulatedData3 <- list(y = data, 
                       lon = lon, lat = lat, locs = locs,
                       a = a, 
                       logb = logb, 
                       logs = logs)

save(simulatedData3, file="simulatedData3.RData")
