library(SpatialExtremes)
set.seed(234)
x <- y <- seq(0, 10, length=20)
coor <- cbind(x, y)
locs <- expand.grid(x,y)
n_loc <- length(x)*length(y)

# Simulate a
## nugget is a constant offset, sill is sigma, range is 2*ell 
a <- rgp(1, coor, cov.mod="powexp",
          mean=60, nugget=0, sill=5, range=10, smooth=1,
          grid=TRUE)
#filled.contour(x, y, a, color.palette = terrain.colors, main="a")

# Simulate b
b <- rgp(1, coor, cov.mod="powexp",
          mean=20, nugget=0, sill=3, range=10, smooth=1,
          grid=TRUE)
logb <- log(b)
#filled.contour(x, y, b, color.palette = terrain.colors, main="b")


# Simulate s 
logs <- rgp(1, coor, cov.mod="powexp",
            mean=-2, nugget=0, sill=1, range=10, smooth=1,
            grid=TRUE)
#filled.contour(x, y, logs, color.palette = terrain.colors, main="s")
s <- exp(logs)

# Simulate data
Y <- Map(rgev, n=sample(10:30, n_loc, replace=TRUE),
         loc=as.vector(a), scale=as.vector(b), shape=as.vector(s))
Y_mat <- matrix(sapply(Y, mean), ncol=sqrt(n_loc))
#filled.contour(x, y, Y_mat, 
#              color.palette = terrain.colors, main="y")

simulatedData2 <- list(locs=locs, a=as.vector(a), logb=as.vector(logb), 
		       logs=logs, y=Y)
save(simulatedData2, file="simulatedData2.RData")
