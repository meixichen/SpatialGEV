set.seed(946)
require(ncdf4) # this package is required to read the original IBTrACS data
require(dplyr)
require(evd)
require(SpatialGEV)

# Import the data and preprocess it
storms <- nc_open("IBTrACS.since1980.v04r00.nc")
lat <- ncvar_get(storms, "lat") 
lon <- ncvar_get(storms, "lon")
time <- ncvar_get(storms, 'iso_time')
mws <- ncvar_get(storms, "wmo_wind") * 0.514 
mws.df <- data.frame(lon=c(lon), lat=c(lat), time=as.Date(c(time)), value=c(mws)) 
mws.df <- mws.df[complete.cases(mws.df), ] # remove NAs
mws.df <- mws.df %>% distinct() %>% arrange(time)
rm(list = c("lat", "lon", "mws", "storms", "time"))

# Explore the distribution of data on a map
map_plot(mws.df$value, mws.df[, c(1,2)])

# Choose a small study region
lon.range <- c(-80, -50)
lat.range <- c(20, 50)
spatial.window <- (mws.df$lon>lon.range[1] & mws.df$lon<=lon.range[2] &
                     mws.df$lat>=lat.range[1] & mws.df$lat<=lat.range[2])
test_set <- mws.df[spatial.window,] # the smaller test set with coordinates in the specified spatial window
map_plot(test_set$value, test_set[, c(1,2)]) # visualize the region

# Grid the data
stations <- grid_location(test_set$lon, test_set$lat, sp.resolution = 2, lon.range = c(-100, -20), lat.range = c(-20, 80))
grid_data <- cbind(test_set, stations)
grid_data <- grid_data %>% group_by(station_ind, station_lon, station_lat) %>% summarise(max_wind = max(value))
map_plot(grid_data$max_wind, as.matrix(grid_data[, c(2,3)]))

# Fit an GEV model to grid_data using MLE. The estimates will be used as initial value in the later optimization of TMB model.
gev_nll <- function(param){
  log_a <- param[1]
  log_b <- param[2]
  log_s <- param[3]
  -sum(dgev(grid_data$max_wind, exp(log_a), exp(log_b), exp(log_s), log = TRUE))
}

test_fit <- optim(par = c(0, 0, 0), gev_nll)
log_param_est <- test_fit$par

# Construct the hierarchical 
y <- grid_data$max_wind
X <- grid_data[,c(2,3)]
theta <- list(a = 50, log_b = rep(log_param_est[2], nrow(X)),
              log_s = log_param_est[3], 
              log_sigma = 0, log_ell = 0)
tmb_model <- hierGev_mod(y = y, x = X, theta = theta, random = "b")

# Optimize the negative log marginal likelihood produced by the model
system.time(
  fit <- optim(tmb_model$par, tmb_model$fn, tmb_model$gr)
)

# Obtain the approximate posterior for fixed and random effects
# p.s. note that the following function took me more than 12 hours to run
draws <- param_sim(tmb_model, y = y, x = X, n_sim = 1000, n_random = nrow(grid_data))
fixed_post_draws <- draws$fixed
random_post_draws <- draws$random

# Plot the approximate posterior
## Fixed effects
par(mfrow = c(2,2))
hist(fixed_post_draws[1, 6:1000], ylab = "", xlab = "a", main = NULL, breaks = 30, 
     cex.axis = 1.2, cex.lab = 1.2) # location
hist(fixed_post_draws[2, 6:1000], ylab = "", xlab = "log(s)", main = NULL, breaks = 30, 
     cex.axis = 1.2, cex.lab = 1.2) # log shape

hist(fixed_post_draws[3, 6:1000], ylab = "", xlab = expression(log(sigma^2)), main = NULL, breaks = 35, 
     cex.axis = 1.2, cex.lab = 1.2) # log sigma
hist(fixed_post_draws[4, 6:1000], ylab = "", xlab= "log(\u2113)", main = NULL, breaks = 30, 
     cex.axis = 1.2, cex.lab = 1.2) # log ell

## Random effect
plot_ind <- sort(sample(1:nrow(grid_data), 9))
par(mfrow = c(3,3))
for (i in plot_ind){
  hist(random_post_draws[i,6:995], breaks = 40,  ylab = "",
       xlab = parse(text = paste0(expression(log(b)), "[", i, "]")),
       main = NULL, cex.axis = 1.6, cex.lab = 2)
}

