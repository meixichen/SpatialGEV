context("model_spde")

random_list <- c("a", "ab", "abs")
reparam_list <- c("positive", "negative", "unconstrained", "zero")

test_that("`model_*_spde` can be compiled on the c++ side.",{
  for (random in random_list){
    for (reparam_s in reparam_list){
      if (random == "abs" & reparam_s == "zero"){
        next
      }
      sim_res <- test_sim(random=strsplit(random, "")[[1]], kernel="spde",
                          reparam_s=reparam_s)
      adfun_spde <- spatialGEV_fit(
        data = sim_res$y,
        locs = sim_res$locs,
        random = random,
        init_param = sim_res$params,
        reparam_s = sim_res$reparam_s,
        kernel = sim_res$kernel,
        adfun_only = TRUE,
        silent = TRUE
      )
    }
  }
  expect_named(adfun_spde, c("adfun", "mesh"))
})
