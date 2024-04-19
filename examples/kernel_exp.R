X1 <- cbind(runif(10, 1, 10), runif(10, 10, 20))
X2 <- cbind(runif(5, 1, 10), runif(5, 10, 20))

kernel_exp(sigma=2, ell=1, X1=X1, X2=X2)

kernel_exp(as.matrix(stats::dist(X1)), sigma=2, ell=1)
