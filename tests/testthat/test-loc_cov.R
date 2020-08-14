context("loc_cov")

test_that("`loc_cov` gives the correct spatial covariance matrix", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    ########### Data simulation #################
    n.obs <- sample(2:10, 1) # nrow of the matrix containing the coordinates
    coor <- cbind(sample(n.obs), sample(n.obs)) # the `n x 2` matrix containing the coordinates
    sig <- exp(rnorm(1))
    ell <- exp(rnorm(1))
    #############################################
    # my covariance matrix 
    mat1 <- matrix(NA, n.obs, n.obs)
    for (i in 1:n.obs){
      for (j in 1:n.obs){
        h <- sqrt(sum((coor[i,]-coor[j,])^2))
        mat1[i,j] <- sig*exp(-h/ell)
      }
    }
    # matrix produced by loc_cov function
    dd <- as.matrix(dist(coor))
    mat2 <- loc_cov(dd, sig, ell)
    # compare the two matrices
    expect_equal(c(mat1), c(mat2))
  }
})