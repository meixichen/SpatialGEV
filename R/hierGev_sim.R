#' Simulate parameters and data to test a hierarchical GEV model.
#' 
#' @param n.obs Number of observations.
#' @return A list of simulated random parameters, a dist object, and response values (MWS).
#' @details Note that the purpose of this package is not to simulate data from the hierarchical model. 
#' Instead, it just conveniently provides a set of data values to test the functions. 
#' @export
hierGev_sim <- function(n.obs){
  # a: location parameter
  # b: scale parameter
  # s: shape parameter
  # sig and ell: hyper parameters for the underlying GP 
  
  # Simulate the hyperparameters for the three underlying GP's
  log.siga <- rnorm(1)
  log.ella <- rnorm(1, mean = 1, sd = 0.2)
  log.sigs <- rnorm(1)
  log.ells <- rnorm(1, mean = 1, sd = 0.2)
  log.sigb <- rnorm(1)
  log.ellb <- rnorm(1, mean = 1, sd = 0.2)
  log.sigs <- rnorm(1)
  log.ells <- rnorm(1, mean = 2, sd = 0.2)
  # Simulate coordinates
  lon <- rnorm(n.obs, mean = 10)
  lat <- rnorm(n.obs, mean = 5)
  x <- cbind(lon,lat)
  colnames(x) <- c("lon","lat")
  loc.mat <- cbind(lon, lat) # an `n x 2` matrix of coordinates
  # calculate the distances between locations
  dd <- as.matrix(dist(loc.mat))
  # simulate the spatial cov matrices
  cov.b <- loc_cov(dd, sigma = exp(log.sigb), ell = exp(log.ellb))
  cov.s <- loc_cov(dd, sigma = exp(log.sigs), ell = exp(log.ells))
  # simulate the scale parameter vector
  log.b <- as.vector(rmvnorm(1, mean = rep(0, n.obs), sigma = cov.b)) 
  # Simulate the shape parameter vector
  log.s <- as.vector(rmvnorm(1, mean = rep(0, n.obs), sigma = cov.s))
  # Simulate the location parameter vector
  y_min <- exp(rnorm(1))
  a.min <- exp(log.s)/exp(log.b) # the lower bound of a
  a.max <- a.min+y_min # the upper bound of a
  a <- mapply(runif, min=a.min+1e-2, max=a.max-1e-2, MoreArgs = list(n=1)) # this way a should guarantee a-b/s>0
  # Simulate response values
  y <- mapply(rgev, loc = a, scale = exp(log.b), shape = exp(log.s), MoreArgs = list(n=1))
  
  # Output
  list(
    a = a,
    log.b = log.b,
    log.s = log.s,
    log.siga = log.siga,
    log.ella = log.ella,
    log.sigb = log.sigb,
    log.ellb = log.ellb,
    log.sigs = log.sigs,
    log.ells = log.ells,
    y_min = y_min,
    dd = dd,
    y = y,
    x = x
  )
}
