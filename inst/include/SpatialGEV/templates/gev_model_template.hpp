#ifndef model_{{random_effects}}_{{kernel}}_hpp
#define model_{{random_effects}}_{{kernel}}_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// TMB specification of GEV-GP models with a chosen covariance kernel.
///
/// The model is defined as follows:
///
/// y ~ GEV(a, b, s),
/// random ~ GP(log_sigma, log_kappa/ell),
///
/// where random can be either one of more of a, b, and s.
///
/// @param[in] y Response vector of length `n_obs`.  Assumed to be > 0.
/// @param[in] loc_ind Location vector of length `n_obs` of integers `0 <= i_loc < n_loc` indicating to which of the SPDE mesh locations each element of `y` is associated.
/// @param[in] reparam_s Integer indicating the type of shape parameter. 0: `s = 0`, i.e., use Gumbel instead of GEV distribution.  1: `s > 0`, in which case we operate on `log(s)`.  2: `s < 0`, in which case we operate on `log(-s)`.  3: unconstrained.
/// @param[in] beta_prior Integer specifying the type of prior on the design matrix coefficients.  1 is weakly informative normal prior and any other numbers means Lebesgue prior `pi(beta) \propto 1`.
/// @param[in] spde Only used for model_x_spde. Object of type `spde_t` as constructed in R by a call to [INLA::inla.spde2.matern()] consisting of `n_loc` mesh locations.

/// @param[in] dd Only used for model_x_exp or model_x_matern. Distance matrix.
/// @param[in] design_mat_x Design matrix of size `n_loc x n_cov`, where `n_cov` is the number of covariates for `x`.
/// @param[in] beta_x_prior Vector of length 2 containing the mean and sd of the normal prior on `beta_x`.
///
/// The following are only used for model_x_matern or matern_x_spde:
/// @param[in] x_pc_prior Integer specifying the type of prior to use on the Matern covariance parameters `log_sigma_a` and `log_kappa_a`.  1 for using PC prior on a, 0 for using Lebesgue prior.
/// @param[in] range_x_prior Vector of length 2 `(rho_0, p_rho)` s.t. `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_x_prior Vector of length 2 `(sig_0, p_sig)` s.t. `Pr(sig > sig_0) = p_sig`.
///
/// The following are only used for model="a" or "ab, i.e., when "s" is treated as a fixed effect:
/// @param[in] s_mean Scalar for Normal prior mean on s.
/// @param[in] s_sd Scalar for Normal prior sd on s.
///
/// @param[in] a GEV location parameter(s).  Vector of length `n_loc` if `is_random_a == 0` or length 1 if `is_random_a == 1`.
/// @param[in] log_b GEV scale parameter(s) on the log scale.  Same shape as `a`.
/// @param[in] s GEV shape parameter(s) on the scale specified by `reparam_s`.  Same shape as `a`.
///
/// @param[in] beta_a, log_sigma_a, log_kappa/ell_a Hyperparameters of Matern GP for `a`.
/// @param[in] beta_b, log_sigma_b, log_kappa/ell_b Hyperparameters of Matern GP for `log_b`.
/// @param[in] beta_s, log_sigma_s, log_kappa/ell_s Hyperparameters of Matern GP for `s`.
template<class Type>
Type model_{{random_effects}}_{{kernel}}(objective_function<Type>* obj){
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;

  // ------ Data inputs ------------
  DATA_VECTOR(y);
  DATA_IVECTOR(loc_ind);
  DATA_INTEGER(reparam_s);
  DATA_INTEGER(beta_prior);
  // SPDE inputs
  {{#use_spde}}
  DATA_STRUCT(spde, spde_t);
  {{/use_spde}}
  {{^use_spde}}
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  {{/use_spde}}
  {{#use_matern}}
  DATA_SCALAR(nu);
  {{/use_matern}}

  // Inputs for a
  {{#is_random_a}}
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);
  {{#use_matern}}
  DATA_INTEGER(a_pc_prior);
  DATA_VECTOR(range_a_prior);
  DATA_VECTOR(sigma_a_prior);
  {{/use_matern}}
  {{/is_random_a}}

  // Inputs for b
  {{#is_random_b}}
  DATA_MATRIX(design_mat_b);
  DATA_VECTOR(beta_b_prior);
  {{#use_matern}}
  DATA_INTEGER(b_pc_prior);
  DATA_VECTOR(range_b_prior);
  DATA_VECTOR(sigma_b_prior);
  {{/use_matern}}
  {{/is_random_b}}

  // Inputs for s
  {{#is_random_s}}
  DATA_MATRIX(design_mat_s);
  DATA_VECTOR(beta_s_prior);
  {{#use_matern}}
  DATA_INTEGER(s_pc_prior);
  DATA_VECTOR(range_s_prior);
  DATA_VECTOR(sigma_s_prior);
  {{/use_matern}}
  {{/is_random_s}}
  {{^is_random_s}}
  DATA_SCALAR(s_mean);
  DATA_SCALAR(s_sd);
  {{/is_random_s}}

  // ------------ Parameters ----------------------
  {{#is_random_a}}
  PARAMETER_VECTOR(a);
  {{/is_random_a}}
  {{^is_random_a}}
  PARAMETER(a);
  {{/is_random_a}}
  {{#is_random_b}}
  PARAMETER_VECTOR(log_b);
  {{/is_random_b}}
  {{^is_random_b}}
  PARAMETER(log_b);
  {{/is_random_b}}
  {{#is_random_s}}
  PARAMETER_VECTOR(s);
  {{/is_random_s}}
  {{^is_random_s}}
  PARAMETER(s);
  {{/is_random_s}}
  {{#is_random_a}}
  PARAMETER_VECTOR(beta_a);
  {{/is_random_a}}
  {{#is_random_b}}
  PARAMETER_VECTOR(beta_b);
  {{/is_random_b}}
  {{#is_random_s}}
  PARAMETER_VECTOR(beta_s);
  {{/is_random_s}}
  {{#is_random_a}}
  PARAMETER({{gp_hyperparam1}}_a);
  PARAMETER({{gp_hyperparam2}}_a);
  {{/is_random_a}}
  {{#is_random_b}}
  PARAMETER({{gp_hyperparam1}}_b);
  PARAMETER({{gp_hyperparam2}}_b);
  {{/is_random_b}}
  {{#is_random_s}}
  PARAMETER({{gp_hyperparam1}}_s);
  PARAMETER({{gp_hyperparam2}}_s);
  {{/is_random_s}}


  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  // ---------- Likelihood contribution from a ------------------
  {{#is_random_a}}
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += nlpdf_gp_{{kernel}}<Type>(mu_a, {{nlpdf_gp_distance}},
                                   exp({{gp_hyperparam1}}_a), exp({{gp_hyperparam2}}_a),
                                   {{nlpdf_gp_extra}});
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior(0), beta_a_prior(1));
  {{#use_matern}}
  nll += nlpdf_matern_hyperpar_prior<Type>({{gp_hyperparam2}}_a, {{gp_hyperparam1}}_a, a_pc_prior,
                                           nu, range_a_prior, sigma_a_prior);
  {{/use_matern}}
  {{/is_random_a}}

  // ---------- Likelihood contribution from b ------------------
  {{#is_random_b}}
  // GP latent layer
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += nlpdf_gp_{{kernel}}<Type>(mu_b, {{nlpdf_gp_distance}},
                                   exp({{gp_hyperparam1}}_b), exp({{gp_hyperparam2}}_b),
                                   {{nlpdf_gp_extra}});
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior(0), beta_b_prior(1));
  {{#use_matern}}
  nll += nlpdf_matern_hyperpar_prior<Type>({{gp_hyperparam2}}_b, {{gp_hyperparam1}}_b, b_pc_prior,
                                           nu, range_b_prior, sigma_b_prior);
  {{/use_matern}}
  {{/is_random_b}}

  // ---------- Likelihood contribution from s ------------------
  {{#is_random_s}}
  // GP latent layer
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll += nlpdf_gp_{{kernel}}<Type>(mu_s, {{nlpdf_gp_distance}},
                                   exp({{gp_hyperparam1}}_s), exp({{gp_hyperparam2}}_s),
                                   {{nlpdf_gp_extra}});
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior(0), beta_s_prior(1));
  {{#use_matern}}
  nll += nlpdf_matern_hyperpar_prior<Type>({{gp_hyperparam2}}_s, {{gp_hyperparam1}}_s, s_pc_prior,
                                           nu, range_s_prior, sigma_s_prior);
  {{/use_matern}}
  {{/is_random_s}}
  {{^is_random_s}}
  // FIXME: rename this to not depend on `s`
  nll += nlpdf_s_prior<Type>(s, s_mean, s_sd);
  {{/is_random_s}}

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), {{a_var}}, {{b_var}}, {{s_var}}, reparam_s);
  }

  {{#calc_z_p}}
  // ------------- Output z -----------------------
  vector<Type> z(loc_ind.size());
  Type p = {{prob}};
  for (int i=0; i<y.size();i++){
    z[i] = {{a_var}}-exp({{b_var}})/{{s_var}}*(1-pow(-log(1-p), -{{s_var}}));
  }
  ADREPORT(z);
  {{/calc_z_p}}

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

