#' Plot a the response values on a map.
#'
#' @param y Vector of `n` response values. 
#' @param x `n x 2`matrix of longitude and latitude of the corresponding response values.
#' @param file A character string naming a file for writing.
#' @param title A character string of the title of the output picture.
#' @param ylim The scale of `y` to generate colors. Default is the range of `y`.
#' @param width The width of the output file. Default is 800.
#' @param height The height of the output file. Default is 600.
#' @return A `png` file of the generated map with different colors indicating the values of y on different locations.
#' @export
map_plot <- function(y, x, file="", title, ylim, width=800, height=600){
  png(file, width=width, height=height, units="px", pointsize=16)
  parset = par(mar=c(0,0,0,5)+2)
  maps::map("world", mar=rep(0,4)+1, xlim=c(-95,-60), ylim=c(15,50), col="white", fill=T, lwd=0.05)
  if (missing(ylim)) ylim <- range(y)
  fields::image.plot(legend.only=T, zlim=ylim, col=viridis(6))
  val = fields::color.scale(y, col=alpha(viridis(6), 1))
  points(as.vector(x[,1]), as.vector(x[,2]), col=val, pch=16, cex=0.6)
  title(title)
  map("world", mar=rep(0,4)+1, xlim=range(x[,1]), ylim=range(x[,2]), lwd=1.25, add=T)
  map("state", mar=rep(0,4)+1, xlim=range(x[,1]),ylim=range(x[,2]), lwd=1.25, add=T)
  par(parset)
  dev.off()
  print("Image generated.")
}
