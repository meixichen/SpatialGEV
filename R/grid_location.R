#' Grid the locations with fixed cell size
#'
#' @param lon Numeric, `n` longitude values
#' @param lat Numeric, `n` latitude values
#' @param sp.resolution Numeric, must be a single value that indicates the minimal unit length of 
#' a grid cell.  
#' @param lon.range Optional vector that indicates the range of `lon`. Default is `range(lon)`.
#' @param lat.range Optional vector that indicates the range of `lat`. Default is `range(lat)`.
#' @return An `n x 3` data frame containing three variables: `cell_ind` corresponds to unique id 
#' for each grid cell, 
#' `cell_lon` is the longitude of the grid cell, `cell_lat` is the latitude of the grid cell. 
#' Since the output data frame retains the order of the input coordinates, the original coordinate 
#' dataset and the output have can be linked one-to-one by the row index.
#' @details The longitude and latitude of each grid cell are the coordinate of the cell center.
#' For example, if `sp.resolution=1`, then `cell_lon=55.5` and `cell_lat=22.5` correspond to the 
#' square whose left boundary is 55, right boundary is 56, upper boundary is 23, and lower boundary
#' is 22.
#' @examples
#' longitude <- runif(20, -90, 80)
#' latitude <- runif(20, 40, 60)
#' grid_locs <- grid_location(longitude, latitude, sp.resolution=0.5)
#' cbind(longitude, latitude, grid_locs)
#' @export

grid_location <- function(lon, lat, sp.resolution=2, lon.range=range(lon), lat.range=range(lat)){
  if (length(lon) != length(lat)){
    stop("lon and lat must have the same length")
  }
  grid.x <- seq(lon.range[1], lon.range[2]-sp.resolution, by = sp.resolution)
  grid.y <- seq(lat.range[1], lat.range[2]-sp.resolution, by = sp.resolution)
  grid.mat <- matrix(1 : (length(grid.x) * length(grid.y)), nrow = length(grid.y))
  x.it <- findInterval(lon, grid.x)
  y.it <- findInterval(lat, grid.y)
  n <- length(lon)
  cell_ind <- rep(0, n)
  cell_lon <- rep(9999, n)
  cell_lat <- rep(9999, n)
  for (i in 1:n){
    cell_ind[i] <- grid.mat[y.it[i], x.it[i]]
    cell_lon[i] <- (grid.x[x.it[i]] + grid.x[x.it[i]] + sp.resolution)/2
    cell_lat[i] <- (grid.y[y.it[i]] + grid.y[y.it[i]] + sp.resolution)/2
  }
  if (any(c(cell_lon==9999, cell_lat==9999))){
    stop("Error occurred in calculating the mean coordinates")
  }
  data.frame(cell_ind, cell_lon, cell_lat)
}
