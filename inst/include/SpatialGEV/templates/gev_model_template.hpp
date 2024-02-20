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
{{#re_names}}
/// {{long_name}} ~ GP({{gp_hyperparam1}}_{{short_name}}, {{gp_hyperparam2}}_{{short_name}})
{{/re_names}}
/// where the GP is parameterized using the {{kernel}} covariance kernel.
///
/// --------- Data provided from R ---------------
/// @param[in] y Response vector of length `n_obs`.  Assumed to be > 0.
/// @param[in] loc_ind Location vector of length `n_obs` of integers `0 <= i_loc < n_loc` indicating
/// to which
/// locations each element of `y` is associated.
/// @param[in] reparam_s Integer indicating the type of shape parameter. 0: `s = 0`, i.e., use
/// Gumbel instead
/// of GEV distribution.  1: `s > 0`, in which case we operate on `log(s)`.  2: `s < 0`, in which
/// case we operate on `log(-s)`.  3: unconstrained.
/// @param[in] beta_prior Integer specifying the type of prior on the design matrix coefficients.
/// 1 is weakly informative normal prior and any other numbers means Lebesgue prior
/// `pi(beta) \propto 1`.
{{#use_spde}}
/// @param[in] spde Object of type `spde_t` as constructed in R by a call to
/// [INLA::inla.spde2.matern()]
/// consisting of `n_loc` mesh locations.
{{/use_spde}}
{{^use_spde}}
/// @param[in] dist_mat `n_loc x n_loc` distance matrix typically constructed via
/// `stats::dist(coordinates)`.
/// @param[in] sp_thres Scalar number used to make the covariance matrix sparse by thresholding.
/// If sp_thres=-1, no thresholding is made.
{{/use_spde}}
{{#re_names}}
/// @param[in] design_mat_{{short_name}} Design matrix of size `n_loc x n_covariate` for parameter
/// {{long_name}}.
/// @param[in] beta_{{short_name}}_prior Vector of length 2 containing the mean and sd of the normal
/// prior on `beta_{{short_name}}`.
{{/re_names}}
{{#use_matern}}
/// @param[in] nu Presepecified smoothness parameter for the Matérn covariance kernel applicable to
/// all random effects.
{{#re_names}}
/// @param[in] {{short_name}}_pc_prior Integer specifying the type of prior to use on the Matérn
/// GP on {{long_name}}.
/// 1 for using PC prior on {{long_name}}, 0 for using Lebesgue prior.
/// `log_sigma_a` and `log_kappa_a`.  1 for using PC prior on a, 0 for using Lebesgue prior.
/// @param[in] range_{{short_name}}_prior PC prior on the range parameter for the Matérn GP on
/// {{long_name}}. Vector of length 2 `(rho_0, p_rho)` s.t. `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_{{short_name}}_prior PC prior on the variance parameter for the Matérn GP on
/// {{long_name}}. Vector of length 2 `(sig_0, p_sig)` s.t. `Pr(sig > sig_0) = p_sig`.
{{/re_names}}
{{/use_matern}}
{{^is_random_s}}
/// @param[in] s_mean Scalar for Normal prior mean on s.
/// @param[in] s_sd Scalar for Normal prior sd on s.
{{/is_random_s}}
///
/// --------- Parameters to estimate ------------
/// @param[in] a GEV location parameter.
{{#is_random_a}}
/// Vector of length `n_loc`.
{{/is_random_a}}
{{^is_random_a}}
/// Vector of length 1.
{{/is_random_a}}
/// @param[in] log_b GEV scale parameter on the log scale.
{{#is_random_b}}
/// Vector of length `n_loc`.
{{/is_random_b}}
{{^is_random_b}}
/// Vector of length 1.
{{/is_random_b}}
/// @param[in] s GEV shape parameter on the scale specified by `reparam_s`.
{{#is_random_s}}
/// Vector of length `n_loc`.
{{/is_random_s}}
{{^is_random_s}}
/// Vector of length 1.
{{/is_random_s}}
{{#re_names}}
/// @param[in] beta_{{short_name}} GP mean covariate coefficient vector of length `n_covariate`
/// for {{long_name}}.
/// @param[in] {{gp_hyperparam1}}_{{short_name}} GP covariance kernel variance hyperparameter
/// for {{long_name}}.
/// @param[in] {{gp_hyperparam2}}_{{short_name}} GP covariance kernel range hyperparameter
/// for {{long_name}}.
{{/re_names}}
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
  {{#use_spde}}
  DATA_STRUCT(spde, spde_t);
  {{/use_spde}}
  {{^use_spde}}
  DATA_MATRIX(dist_mat);
  DATA_SCALAR(sp_thres);
  {{/use_spde}}
  {{#use_matern}}
  DATA_SCALAR(nu);
  {{/use_matern}}

  {{#re_names}}
  // Inputs for {{long_name}}
  DATA_MATRIX(design_mat_{{short_name}});
  DATA_VECTOR(beta_{{short_name}}_prior);
  {{#use_matern}}
  DATA_INTEGER({{short_name}}_pc_prior);
  DATA_VECTOR(range_{{short_name}}_prior);
  DATA_VECTOR(sigma_{{short_name}}_prior);
  {{/use_matern}}
  {{/re_names}}
  {{^is_random_s}}
  DATA_SCALAR(s_mean);
  DATA_SCALAR(s_sd);
  {{/is_random_s}}

  // ------------ Parameters ----------------------

  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);

  {{#re_names}}
  PARAMETER_VECTOR(beta_{{short_name}});
  {{/re_names}}
  {{#re_names}}
  PARAMETER({{gp_hyperparam1}}_{{short_name}});
  PARAMETER({{gp_hyperparam2}}_{{short_name}});
  {{/re_names}}

  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  {{#re_names}}
  // ---------- Likelihood contribution from {{long_name}} ------------------
  // GP latent layer
  vector<Type> mu_{{short_name}} = {{long_name}} - design_mat_{{short_name}} * beta_{{short_name}};
  nll += nlpdf_gp_{{kernel}}<Type>(mu_{{short_name}}, {{nlpdf_gp_distance}},
				   exp({{gp_hyperparam1}}_{{short_name}}),
				   exp({{gp_hyperparam2}}_{{short_name}}),
                                   {{nlpdf_gp_extra}});
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_{{short_name}}, beta_prior, beta_{{short_name}}_prior(0),
                                beta_{{short_name}}_prior(1));
  {{#use_matern}}
  nll += nlpdf_matern_hyperpar_prior<Type>({{gp_hyperparam2}}_{{short_name}},
					   {{gp_hyperparam1}}_{{short_name}},
					   {{short_name}}_pc_prior,
                                           nu, range_{{short_name}}_prior,
					   sigma_{{short_name}}_prior);
  {{/use_matern}}
  {{/re_names}}
  {{^is_random_s}}
  // FIXME: rename this to not depend on `s`
  nll += nlpdf_s_prior<Type>({{s_var_loc}}, s_mean, s_sd);
  {{/is_random_s}}

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), {{a_var_loc}}, {{b_var_loc}}, {{s_var_loc}}, reparam_s);
  }

  {{#calc_z_p}}
  // ------------- Output z -----------------------
  DATA_INTEGER(return_level);
  vector<Type> z({{re}}.size());
  if (return_level == 1){
    Type p = {{prob}};
    for (int i=0; i<{{re}}.size();i++){
      z[i] = {{a_var}}-exp({{b_var}})/{{s_var}}*(1-pow(-log(1-p), -{{s_var}}));
    }
  }
  ADREPORT(z);
  {{/calc_z_p}}

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif


