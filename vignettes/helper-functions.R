#' Plot a quantity on a grid map.
#'
#' @param x Vector of longitude.
#' @param y Vector of latitude.
#' @param z Matrix of the quantity to plot.
#' @param title Character string.
#' @param x_lab Character string for x-axis label.
#' @param y_lab Character string for y-axis label.
#' @param cex Font size.
#'
#' @return A plot object.
grid_plot <- function(x, y, z,
                      title, x_lab="Longitude", y_lab="Latitude", cex=1.2) {
  fields::image.plot(x=x, y=y, z=z,
             xlab=x_lab, ylab=y_lab, main=title,
             cex.lab=cex, cex.axis=cex, axis.args=list(cex.axis=cex))
}

#' Plot values on a map.
#'
#' @param value Vector of values to plot.
#' @param zlim Range of the values to plot.
#' @param lon Vector of longitudes.
#' @param lat Vector of latitudes.
#' @param title Character string.
#'
#' @return A plot object.
map_plot <- function(value, zlim=NULL,
                     lon=locs$cell_lon, lat=locs$cell_lat, title="") {
  if (is.null(zlim)){
    zlim <- range(value)
  }
  val <- fields::color.scale(value, col=viridisLite::viridis(10),
                             zlim=zlim)
  plot(lon, lat, col=val, axes=FALSE, pch=15,
       cex=1.2, xlab="", ylab="", main=title)
  maps::map("world", "Canada", add=TRUE)
  fields::image.plot(legend.only = TRUE,
                     zlim=zlim, col=viridisLite::viridis(10),
                     horizontal = FALSE, legend.shrink = 0.5,
                     cex=0.7, smallplot=c(0.85, 0.9, 0.5, 0.9))
  axis(1, at=seq(-140, -50, by=10))
  axis(2, at=seq(40, 80, by=10))
}

