#######################################################
# Import and preprocess the IBTrACS dataset
# Data available at https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/netcdf/
# Knapp, K. R., M. C. Kruk, D. H. Levinson, H. J. Diamond, and C. J. Neumann, 2010: The International Best Track Archive for Climate Stewardship (IBTrACS): Unifying tropical cyclone best track data
######################################################
require(ncdf4) # this package is required to read the original IBTrACS data
require(dplyr)
require(SpatialGEV)
storms <- nc_open("IBTrACS.since1980.v04r00.nc")

lat <- ncvar_get(storms, "lat") 
lon <- ncvar_get(storms, "lon")
time <- ncvar_get(storms, 'iso_time')
mws <- ncvar_get(storms, "wmo_wind") * 0.514 
mws.df <- data.frame(lon=c(lon), lat=c(lat), time=as.Date(c(time)), value=c(mws)) 
mws.df <- mws.df[complete.cases(mws.df), ] # remove NAs
mws.df <- mws.df %>% distinct() %>% arrange(time)
rm(list = c("lat", "lon", "mws", "storms", "time"))
lon.range <- c(-120, 0)
lat.range <- c(10, 40)
spatial.window <- (mws.df$lon>lon.range[1] & mws.df$lon<=lon.range[2] &
                     mws.df$lat>=lat.range[1] & mws.df$lat<=lat.range[2])
test_set <- mws.df[spatial.window,] # the smaller test set with coordinates in the specified spatial window

# Grid the data
stations <- grid_location(test_set$lon, test_set$lat, sp.resolution = 3, lon.range = lon.range, lat.range = lat.range)
grid_data <- cbind(test_set, stations)
grid_data <- grid_data %>% group_by(station_ind, station_lon, station_lat) %>% summarise(max_wind = max(value), count=n())
nrow(grid_data)
grid_data <- grid_data[which(grid_data$count>=20),]
X <- grid_data %>% ungroup() %>% select(station_lon, station_lat)
dd.mat <- as.matrix(dist(X))
n_obs <- nrow(grid_data)

## Plot
# pdf("case-study-raw-map.pdf", width=15, height=8)
# map_plot(test_set$value, test_set[, c(1,2)], legend.shrink = 0.5, legend.horizontal = TRUE, legend.position = c(0.52, 0.9, 0.2,0.22)) # visualize the region
# dev.off()
# pdf("case-study-grid-map.pdf", width=15, height=9)
# map_plot(grid_data$max_wind, as.matrix(grid_data[, c(2,3)]), legend.horizontal = TRUE, legend.position = c(0.6, 0.95, 0.2,0.22))