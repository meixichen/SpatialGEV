require(whisker)
require(TMBtools)
require(usethis)
pkg_dir <- usethis::proj_get()
template_file <- file.path(pkg_dir, 
                           "inst", "include", "SpatialGEV", "templates", 
                           "gev_model_template.hpp")
template <- readLines(template_file)

#---------- Helper functions for parsing the template ----------------
choose_gp_hyperparam <- function(kernel = c("exp", "matern", "spde")){
  kernel <- match.arg(kernel)
  switch(kernel,
         exp = c("log_sigma", "log_ell"),
         matern = c("log_sigma", "log_kappa"),
         spde = c("log_sigma", "log_kappa"))
}
choose_abs_var_name <- function(random_effects = c("a", "ab", "abs"),
                                with_loc_ind = F){
  random_effects <- match.arg(random_effects)
  out <- switch(random_effects, 
                a = c("a(i)", "log_b(0)", "s(0)"),
                ab = c("a(i)", "log_b(i)", "s(0)"),
                abs = c("a(i)", "log_b(i)", "s(i)"))
  out
}
choose_nlpdf_gp_setting <- function(kernel = c("exp", "matern", "spde")){
  kernel <- match.arg(kernel)
  switch(kernel,
         exp = c("dist_mat", "sp_thres"),
         matern = c("dist_mat", "nu, sp_thres"),
         spde = c("spde", "nu"))
}
create_re_long_short_names <- function(re_logical = c(TRUE, TRUE, TRUE)){
  out <- list(c(short_name="a", long_name="a"),
              c(short_name="b", long_name="log_b"),
              c(short_name="s", long_name="s"))
  out[re_logical] 
}

# ------------- Generate all model combinations -------------------------
# Specify the model to use
random_effects_list <- c("a", "ab", "abs")
kernel_list <- c("exp", "matern", "spde")
re_kernel_combs <- expand.grid(random_effects_list, kernel_list,
                               stringsAsFactors = F)
colnames(re_kernel_combs) <- c("random", "kernel")

for (i in 1:nrow(re_kernel_combs)){
  random_effects <- re_kernel_combs[i, "random"]
  kernel <- re_kernel_combs[i, "kernel"]
  # Keys for generating the template
  check_random_abs <- unname(parse_random(random_effects))
  gp_hyperparam <- choose_gp_hyperparam(kernel)
  abs_var_name <- choose_abs_var_name(random_effects)
  abs_var_name_loc <- gsub("i", "loc_ind(i)", abs_var_name)
  nlpdf_gp_setting <- choose_nlpdf_gp_setting(kernel)
  re_names <- create_re_long_short_names(check_random_abs)
  temp_keys <- list(
    re_names = re_names,
    is_random_a = check_random_abs[1],
    is_random_b = check_random_abs[2],
    is_random_s = check_random_abs[3],
    random_effects = random_effects,
    kernel = match.arg(kernel, c("exp", "matern", "spde")),
    use_spde = kernel=="spde",
    use_matern = kernel %in% c("matern", "spde"),
    a_var_loc = abs_var_name_loc[1],
    b_var_loc = abs_var_name_loc[2],
    s_var_loc = abs_var_name_loc[3],
    nlpdf_gp_distance = nlpdf_gp_setting[1],
    nlpdf_gp_extra = nlpdf_gp_setting[2],
    gp_hyperparam1 = gp_hyperparam[1],
    gp_hyperparam2 = gp_hyperparam[2],
    #calc_z_p = F
    calc_z_p = list(re = unname(re_names[[1]]['long_name']),
                    a_var = abs_var_name[1],
                    b_var = abs_var_name[2],
                    s_var = abs_var_name[3],
                    prob = 0.1)
  )

  write_dir <- file.path(pkg_dir, "src", "TMB")
  writeLines(whisker.render(template, temp_keys),
             file.path(write_dir,
                       paste0(paste("model", random_effects, kernel, sep = "_"), 
                              ".hpp")))
}



TMBtools::export_models()
