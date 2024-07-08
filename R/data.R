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


#' Simulated dataset 1 with covariates
#'
#' A list of data used for package testing and demos. Both `a` and `logb` are simulated on smooth
#' deterministic surfaces depending linearly on two covariates.
#' @format A list containing the simulation parameters and simulated data on a 20x20 grid:
#' \describe{
#'   \item{locs}{A 400x2 matrix. First column contains longitudes and second contains latitudes}
#'   \item{a}{A length 400 vector. GEV location parameters}
#'   \item{logb}{A length 400 vector. Log-transformed GEV scale parameters}
#'   \item{logs}{A scalar. Log-transformed GEV shape parameter shared across space}
#'   \item{y}{A length 400 list of vectors which are observations simulated at
#'   each location}
#'   \item{beta_a}{A length 3 vector containing the regression coefficients for
#'   a, with the first element being the intercept, the second element as the
#'   coefficient for covariate 1, and the third for covariate 2}
#'   \item{beta_b}{A length 2 vector containing the regression coefficients for
#'   b, with the first element being the intercept and the second element as the
#'   coefficient for covariate 1}
#'   \item{covariates}{A 400x3 design matrix containing the covariates, with the
#'   first column being all 1s for intercept}
#' }
"simulatedData_withcov"


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

#' Gridded monthly total snowfall in Canada from 1987 to 2021.
#'
#' Variables containing the monthly total snowfall (in cm) in Canada from 1987 to
#' 2021 and the location information. The data has been gridded and information about
#' the grid size can be found in the paper Fast and Scalable Inference for Spatial
#' Extreme Value Models (arxiv: 2110.07051).
#'
#' @format A list containing the location information and the observations:
#' \describe{
#'   \item{locs}{A 509x2 matrix with longitude and latitude for each grid cell}
#'   \item{n_loc}{Number of locations}
#'   \item{Y}{A list of length 509 with each element of the list containing the
#'   observations at a location}
#' }
#' @source \url{https://climate-change.canada.ca/climate-data/#/monthly-climate-summaries}
"CAsnow"
