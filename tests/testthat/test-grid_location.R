context("grid_location")

test_that("`grid_location` produces the correct grid coordinates.",{
  for (ii in seq_len(50)){
    n_loc <- sample(2:20, 1)
    longitude <- runif(n_loc, runif(1, 1, 9), runif(1, 11, 20))
    latitude <- runif(n_loc, runif(1, 1, 9), runif(1, 11, 20))
    sp_resl <- runif(1, 0, 1)
    grid_locs <- grid_location(longitude, latitude, sp.resolution=sp_resl)
    # Lon and lat of the location should be within sqrt(sp_resl^2/2) distance
    # from lon and lat of the cell center
    expect_true(max(abs(grid_locs$cell_lon-longitude)) < sqrt(sp_resl^2/2))
    expect_true(max(abs(grid_locs$cell_lat-latitude)) < sqrt(sp_resl^2/2))
  }
})

test_that("`grid_location` gives the same cell id to coordinates in the same
          grid cell.", {
  for (ii in seq_len(50)){
    n_loc <- sample(2:20, 1)
    sp_resl <- runif(1, 0, 1)
    lon_start <- runif(1, 0, 10)
    longitude <- runif(n_loc, lon_start, runif(1, lon_start, lon_start+sp_resl))
    latitude <- runif(n_loc, runif(1, 1, 10), runif(1, 11, 20))
    grid_locs <- grid_location(longitude, latitude, sp.resolution=sp_resl)
    expect_equal(length(unique(grid_locs$cell_lon)), 1)

    latitude <- runif(n_loc, lon_start, runif(1, lon_start, lon_start+sp_resl))
    grid_locs <- grid_location(longitude, latitude, sp.resolution=sp_resl)
    expect_equal(unique(grid_locs), grid_locs[1,])
  }
})

test_that("`grid_location` gives errors for incorrect input dimensions.",{
  for (ii in seq_len(50)){
    expect_error(grid_location(runif(1, 1, 10), runif(1, 1, 10),
                               sp.resolution=runif(1, 0, 1)),
                 "There should be at least 2 locations.")

    n_locs <- sample(2:40, 2, replace = FALSE)
    longitude <- runif(n_locs[1], 1, 10)
    latitude <- runif(n_locs[2], 1, 10)
    expect_error(grid_location(longitude, latitude,
                               sp.resolution=runif(1, 0, 1)),
                 "lon and lat must have the same length.")
  }
})
