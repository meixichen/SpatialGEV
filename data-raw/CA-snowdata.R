require(dplyr)
fnames <- list.files(pattern="climate-monthly-*")
i <- 1
fname <- fnames[i]
allsnowdata <- read.csv(fname) %>%
  select(LATITUDE, LONGITUDE, STATION_NAME, CLIMATE_IDENTIFIER,
         PROVINCE_CODE, LOCAL_YEAR, LOCAL_MONTH, TOTAL_SNOWFALL) %>%
  na.omit()

for (i in 2:length(fnames)){
  fname <- fnames[i]
  allsnowdata <- read.csv(fname) %>%
    select(LATITUDE, LONGITUDE, STATION_NAME, CLIMATE_IDENTIFIER,
           PROVINCE_CODE, LOCAL_YEAR, LOCAL_MONTH, TOTAL_SNOWFALL) %>%
    na.omit() %>%
    full_join(allsnowdata)
}

raw_locs <- allsnowdata %>%
  filter(TOTAL_SNOWFALL > 0) %>%
  select(LATITUDE, LONGITUDE) %>%
  distinct()
pdf("case-study-data-map-raw.pdf", width=6, height=4.5)
par(mar=c(0.5,1,0.5,1))
plot(raw_locs$LONGITUDE[-which.min(raw_locs[,1])], raw_locs$LATITUDE[-which.min(raw_locs[,1])],
     axes=F, pch=20, cex=0.3, col="red")
maps::map("world", "Canada", add=T)
dev.off()

# Gridding
grid_size <- 1
grid_locs <- grid_location(allsnowdata$LONGITUDE, allsnowdata$LATITUDE,
                           sp.resolution = grid_size)
data_grid <- cbind(grid_locs, allsnowdata)
# Yearly max for each location
all_locs <- data_grid %>%
  select(cell_ind, cell_lon, cell_lat) %>%
  distinct()
yearly_max_records <- data_grid %>%
  group_by(cell_ind, LOCAL_YEAR) %>%
  slice(YEARLY_MAX_SNOWFALL = which.max(TOTAL_SNOWFALL)) %>%
  select(cell_ind, LOCAL_YEAR, LOCAL_MONTH, TOTAL_SNOWFALL) %>%
  rename(YEARLY_MAX_SNOWFALL = TOTAL_SNOWFALL) %>%
  filter(YEARLY_MAX_SNOWFALL > 0) %>% # Remove records of 0s 
  left_join(all_locs, by="cell_ind")

# Coordinates of the locations
locs <- yearly_max_records %>% ungroup() %>%
  select(cell_ind, cell_lon, cell_lat) %>%
  distinct()
n_loc <- nrow(locs)

# Make data into a list in which each vector contains data from one location
Y <- vector(mode="list", length=n_loc)
for (i in 1:n_loc){
  id <- locs$cell_ind[i]
  Y[[i]] <- yearly_max_records %>%
    ungroup() %>%
    filter(cell_ind==id) %>%
    pull(YEARLY_MAX_SNOWFALL)
}

# Only keep locations with at least T years of records
T <- 10
chosen_loc_ind <- which(sapply(Y, length) >= T)
Y <- Y[chosen_loc_ind]
locs <- locs %>% select(cell_lon, cell_lat) %>% slice(chosen_loc_ind)
n_loc <- nrow(locs)
save(Y, locs, n_loc, file="CA-snowdata.RData")

