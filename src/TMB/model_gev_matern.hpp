#ifndef SpatialGEV_model_gev_matern_hpp
#define SpatialGEV_model_gev_matern_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// TMB specification of GEV-GP models with Matern covariance.
///
/// The model is defined as follows:
///
/// [TBD]
///
///
/// @param[in] y Response vector of length `n_obs`.  Assumed to be > 0.
/// @param[in] loc_ind Location vector of length `n_obs` of integers `0 <= i_loc < n_loc` to which each observation in `y` is associated.
///
/// @param[in] reparam_s Integer indicating the type of shape parameter. 0: `s = 0`, i.e., use Gumbel instead of GEV distribution.  1: `s > 0`, in which case we operate on `log(s)`.  2: `s < 0`, in which case we operate on `log(-s)`.  3: unconstrained.
///
/// @param[in] beta_prior Integer specifying the type of prior on the design matrix coefficients.  1 is weakly informative normal prior and any other numbers means Lebesgue prior `pi(beta) \propto 1`.
///
/// @param[in] is_random_a Integer 0 or 1 as to whether parameter a is modeled as a random or fixed effect.  DATA inputs `design_mat_a`, `beta_a_prior`, `a_pc_prior`, `range_a_prior`, `sigma_a_prior` and PARAMETER inputs `beta_a`, `log_sigma_a`, and `log_kappa_a` are unused if `is_random_a == true`.  DATA input `a_prior` is unsued if `is_random_a = false`.
///
/// @param[in] design_mat_a Design matrix of size `n_loc x n_cov_a`.
/// @param[in] beta_a_prior Vector of length 2 containing the mean and sd of the normal prior on `beta_a`.
///
/// @param[in] a_pc_prior Integer specifying the type of prior to use on the Matern covariance parameters `log_sigma_a` and `log_kappa_a`.  1 for using PC prior on a, 0 for using Lebesgue prior.
/// @param[in] range_a_prior Vector of length 2 `(rho_0, p_rho)` s.t. `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_a_prior Vector of length 2 `(sig_0, p_sig)` s.t. `Pr(sig > sig_0) = p_sig`.
/// @param[in] a_fixed_prior Vector of length 2 giving the prior mean and standard deviation of `a` if it is a fixed parameter, i.e., if `is_random_a == false`.  If `a_fixed_prior(1) > 9999` use a Lebesgue prior.
///
/// @param[in] is_random_b, ..., b_fixed_prior Same as for `*_a` but for `log_b`.
/// @param[in] is_random_s, ..., s_fixed_prior Same as for `*_a` but for `log_s`.
///
/// @param[in] use_spde An integer indicating whether to use the SPDE approximation.
/// @param[in] nu Smoothness parameter for Matern covariance (same for all random effects).
/// @param[in] spde Optional object of type `spde_t` as constructed in R by a call to [INLA::inla.spde2.matern()].  Unused if `use_spde == false`.
///
/// @param[in] a GEV location parameter(s).  Vector of length `n_loc` if `is_random_a == 0` or length 1 if `is_random_a == 1`.
/// @param[in] log_b GEV scale parameter(s) on the log scale.  Same shape as `a`.
/// @param[in] s GEV shape parameter(s) on the scale specified by `reparam_s`.  Same shape as `a`.
///
/// @param[in] beta_a, log_sigma_a, log_kappa_a Hyperparameters of Matern GP for `a`.
/// @param[in] beta_b, log_sigma_b, log_kappa_b Hyperparameters of Matern GP for `log_b`.
/// @param[in] beta_s, log_sigma_s, log_kappa_s Hyperparameters of Matern GP for `s`.
template<class Type>
Type model_gev_matern(objective_function<Type>* obj) {
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y);
  DATA_IVECTOR(loc_ind);
  DATA_INTEGER(reparam_s); 
  DATA_INTEGER(beta_prior); 

  // spde inputs
  DATA_INTEGER(use_spde);
  DATA_SCALAR(nu); 
  DATA_STRUCT(spde, spde_t); 

  // a inputs
  DATA_INTEGER(is_random_a);
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);
  DATA_INTEGER(a_pc_prior); 
  DATA_VECTOR(range_a_prior);
  DATA_VECTOR(sigma_a_prior);
  DATA_VECTOR(a_fixed_prior);

  // b inputs
  DATA_INTEGER(is_random_b);
  DATA_MATRIX(design_mat_b);
  DATA_VECTOR(beta_b_prior);
  DATA_INTEGER(b_pc_prior); 
  DATA_VECTOR(range_b_prior);
  DATA_VECTOR(sigma_b_prior);
  DATA_VECTOR(b_fixed_prior);

  // s inputs
  DATA_INTEGER(is_random_s);
  DATA_MATRIX(design_mat_s);
  DATA_VECTOR(beta_s_prior);
  DATA_INTEGER(s_pc_prior); 
  DATA_VECTOR(range_s_prior);
  DATA_VECTOR(sigma_s_prior);
  DATA_VECTOR(s_fixed_prior);
  
  // parameter list
  PARAMETER_VECTOR(a);  
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);

  // hyperparameters for a
  PARAMETER_VECTOR(beta_a);
  PARAMETER(log_sigma_a);
  PARAMETER(log_kappa_a);

  // hyperparameters for b
  PARAMETER_VECTOR(beta_b);
  PARAMETER(log_sigma_b);
  PARAMETER(log_kappa_b);

  // hyperparameters for s
  PARAMETER_VECTOR(beta_s);
  PARAMETER(log_sigma_s);
  PARAMETER(log_kappa_s);
  
   // calculate the negative log likelihood
  Type nll = Type(0.0);

  // data layer
  Type _a, _log_b, _s;
  for(int i=0;i<y.size();i++) {
    _a = a(is_random_a ? loc_ind(i) : 0);
    _log_b = log_b(is_random_b ? loc_ind(i) : 0);
    _s = s(is_random_s ? loc_ind(i) : 0);
    nll -= gev_reparam_lpdf<Type>(y(i), _a, _log_b, _s, reparam_s);
  }

  // GP latent layer
  if(is_random_a) {
    vector<Type> mu_a = a - design_mat_a * beta_a;
    nll += nlpdf_gp_spde<Type>(mu_a, spde, exp(log_sigma_a), exp(log_kappa_a), nu);
  }
  if(is_random_b) {
    vector<Type> mu_b = log_b - design_mat_b * beta_b;
    nll += nlpdf_gp_spde<Type>(mu_b, spde, exp(log_sigma_b), exp(log_kappa_b), nu);
  }
  if(is_random_s) {
    vector<Type> mu_s = s - design_mat_s * beta_s;
    nll += nlpdf_gp_spde<Type>(mu_s, spde, exp(log_sigma_s), exp(log_kappa_s), nu);
  }

  // prior
  if(is_random_a) {
    nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior(0), beta_a_prior(1));
    nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_a, log_sigma_a, a_pc_prior,
					     nu, range_a_prior, sigma_a_prior);
  } else {
    // FIXME: rename this to not depend on `s`
    nll += nlpdf_s_prior<Type>(a(0), a_fixed_prior(0), a_fixed_prior(1));
  }
  if(is_random_b) {
    nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior(0), beta_b_prior(1));
    nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_b, log_sigma_b, b_pc_prior,
					     nu, range_b_prior, sigma_b_prior);
  } else {
    nll += nlpdf_s_prior<Type>(log_b(0), b_fixed_prior(0), b_fixed_prior(1));
  }
  if(is_random_s) {
    nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior(0), beta_s_prior(1));
    nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_s, log_sigma_s, s_pc_prior,
					     nu, range_s_prior, sigma_s_prior);
  } else {
    nll += nlpdf_s_prior<Type>(s(0), s_fixed_prior(0), s_fixed_prior(1));
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

