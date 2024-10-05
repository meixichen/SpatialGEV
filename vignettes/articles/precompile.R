knitr::knit(input = "_SpatialGEV.Rmd", output = "SpatialGEV.Rmd")
file.copy(
  from = "SpatialGEV.Rmd",
  to = "../SpatialGEV.Rmd",
  overwrite = TRUE
)
file.copy(
  from = "SpatialGEV_figures",
  to = "../",
  overwrite = TRUE,
  recursive = TRUE
)
rmarkdown::render("../SpatialGEV.Rmd")
