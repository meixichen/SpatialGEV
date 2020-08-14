context("map_plot")

test_that("`map_plot` is able to generate pdf image output.", {
  n <- 100
  y <- rnorm(n)
  x <- cbind(runif(n, -100,-60), runif(n, 15, 50))
  dir <- paste(getwd(), "test.png", sep = "/")
  expect_output(map_plot(y=y, x=x, file = dir, title = "test"), "Image generated.")
})  