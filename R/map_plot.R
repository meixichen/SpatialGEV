#' Plot values with coordinates on a map.
#'
#' @param y Vector of `n` response values. 
#' @param x `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param title Optional. A character string of the title of the output picture.
#' @param ylim The scale of `y` to generate colors. Default is the range of `y`.
#' @param legend.shrink Relative size of the legend. Default is 0.5.
#' @param legend.horizontal Whether the legend is horizontal in the plot. Default is `FALSE`.
#' @param legend.position Vector of 4 coordinates of the legend each ranging from 0 to 1
#' @return A plot of the generated map with different colors indicating the values of y on different locations.
#' @export
map_plot <- function(y, x, title, ylim, legend.shrink=0.5, legend.horizontal=FALSE, legend.position=NULL){
  if (missing(ylim)) ylim <- range(y)
  maps::map("world", mar=rep(0,4), xlim=range(x[,1]), ylim=range(x[,2]), col="white", fill=T, lwd=1.25)
  maps::map("state", mar=rep(0,4), xlim=range(x[,1]),ylim=range(x[,2]), lwd=1.25, add=T)
  fields::image.plot(legend.only=T, zlim=ylim, col=viridisLite::viridis(6), horizontal = legend.horizontal,
                     legend.shrink = legend.shrink, cex = 1.4, smallplot = legend.position)
  val = fields::color.scale(y, col=scales::alpha(viridisLite::viridis(6), 1))
  points(as.vector(x[,1]), as.vector(x[,2]), col=val, pch=16, cex=0.6)
  if (!missing(title)) title(title)
}
