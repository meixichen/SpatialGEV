longitude <- runif(20, -90, 80)
latitude <- runif(20, 40, 60)
grid_locs <- grid_location(longitude, latitude, sp.resolution=0.5)
cbind(longitude, latitude, grid_locs)
