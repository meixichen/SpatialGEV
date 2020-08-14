#' Grid the locations.
#'
#' @param lon Numeric, `n` longitude values
#' @param lat Numeric, `n` latitude values
#' @param sp.resolution Numeric, must be a single value that indicates the minimal unit length of the grid.  
#' @param lon.range Optional vector that indicates the range of `lon`. Default is `range(lon)`.
#' @param lat.range Optional vector that indicates the range of `lat`. Default is `range(lat)`.
#' @return An `n x 3` data frame containing three variables: `station_ind` corresponds to unique numbers for each grid, `station_lon` is the longitude of the grid, `station_lat` is the latitude of the grid.
#' @details The longitude and latitude of each grid is calculated by averaging the lon/lat of all locations in the same grid.
#' @export

grid_location <- function(lon, lat, sp.resolution=2, lon.range=range(lon), lat.range=range(lat)){
  #require(dplyr)
  if (length(lon) != length(lat)){
    stop("lon and lat must have the same length")
  }
  grid.x <- seq(lon.range[1], lon.range[2]-sp.resolution, by = sp.resolution)
  grid.y <- seq(lat.range[1], lat.range[2]-sp.resolution, by = sp.resolution)
  grid.mat <- matrix(1 : (length(grid.x) * length(grid.y)), nrow = length(grid.y))
  rownames(grid.mat) <- grid.y
  colnames(grid.mat) <- grid.x
  x.it <- findInterval(lon, grid.x)
  y.it <- findInterval(lat, grid.y)
  n <- length(lon)
  station_ind <- rep(0, n)
  for (i in 1:n){
    station_ind[i] <- grid.mat[y.it[i], x.it[i]]
  }
  temp <- data.frame(lon, lat, station_ind)
  station_coor <- temp %>% group_by(station_ind) %>% summarise(lon_mean = mean(lon), lat_mean = mean(lat))
  station_lon <- rep(9999, n)
  station_lat <- rep(9999, n)
  for (i in 1:n){
    ind <- which(station_coor$station_ind == station_ind[i])
    station_lon[i] <- station_coor$lon_mean[ind]
    station_lat[i] <- station_coor$lat_mean[ind]
  }
  if (any(c(station_lon==9999, station_lat==9999))){
    stop("Error occurred in calculating the mean coordinates")
  }
  data.frame(station_ind, station_lon, station_lat)
}