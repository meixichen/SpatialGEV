#' Monthly total snowfall in Ontario, Canada from 1987 to 2021.
#'
#' A dataset containing the monthly total snowfall (in cm) in Ontario, Canada from 1987 to 2021.
#'
#' @format A data frame with 63945 rows and 7 variables with each row corresponding to a monthly
#' record at a weather location:
#' \describe{
#'   \item{LATITUDE}{Numeric. Latitude of the weather station}
#'   \item{LONGITUDE}{Numeric. Longitude of the weather station}
#'   \item{STATION_NAME}{Character. Name of the weather station}
#'   \item{CLIMATE_IDENTIFIER}{Character. Unique id of each station}
#'   \item{LOCAL_YEAR}{Integer from 1987 to 2021. Year of the record}
#'   \item{LOCAL_MONTH}{Integer from 1 to 12. Month of the record}
#'   \item{TOTAL_SNOWFALL}{Positive number. Total monthly snowfall at a station in cm}
#' }
#' @source \url{https://climate-change.canada.ca/climate-data/#/monthly-climate-summaries}
"ONsnow"


#' Simulated dataset 1
#'
#' A list of data used for package testing and demos. Both `a` and `logb` are simulated on smooth
#' deterministic surfaces.
#' @format A list containing the simulation parameters and simulated data on a 20x20 grid:
#' \describe{
#'   \item{locs}{A 400x2 matrix. First column contains longitudes and second contains latitudes}
#'   \item{a}{A length 400 vector. GEV location parameters}
#'   \item{logb}{A length 400 vector. Log-transformed GEV scale parameters}
#'   \item{logs}{A scalar. Log-transformed GEV shape parameter shared across space}
#'   \item{y}{A length 400 list of vectors which are observations simulated at each location}
#' }
"simulatedData"


#' Simulated dataset 2
#'
#' A list of data used for package testing and demos. `a`, `logb`, `logs` are simulated from 
#' respective Gaussian random fields and thus are nonsmooth.
#' @format A list containing the simulation parameters and simulated data on a 20x20 grid:
#' \describe{
#'   \item{locs}{A 400x2 matrix. First column contains longitudes and second contains latitudes}
#'   \item{a}{A length 400 vector. GEV location parameters}
#'   \item{logb}{A length 400 vector. Log-transformed GEV scale parameters}
#'   \item{logs}{A length 400 vector. Log-transformed GEV shape parameters}
#'   \item{y}{A length 400 list of vectors which are observations simulated at each location}
#' }
"simulatedData2"
