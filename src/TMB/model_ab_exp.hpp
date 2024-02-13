#ifndef model_ab_exp_hpp
#define model_ab_exp_hpp

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
Type model_ab_exp(objective_function<Type>* obj){
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
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.

  // Inputs for a
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);

  // Inputs for b
  DATA_MATRIX(design_mat_b);
  DATA_VECTOR(beta_b_prior);

  // Inputs for s
  DATA_SCALAR(s_mean);
  DATA_SCALAR(s_sd);

  // ------------ Parameters ----------------------
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER(s);
  PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(beta_b);
  PARAMETER(log_sigma_a);
  PARAMETER(log_ell_a);
  PARAMETER(log_sigma_b);
  PARAMETER(log_ell_b);


  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  // ---------- Likelihood contribution from a ------------------
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += nlpdf_gp_exp<Type>(mu_a, dd,
                                   exp(log_sigma_a), exp(log_ell_a),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior(0), beta_a_prior(1));

  // ---------- Likelihood contribution from b ------------------
  // GP latent layer
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += nlpdf_gp_exp<Type>(mu_b, dd,
                                   exp(log_sigma_b), exp(log_ell_b),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior(0), beta_b_prior(1));

  // ---------- Likelihood contribution from s ------------------
  // FIXME: rename this to not depend on `s`
  nll += nlpdf_s_prior<Type>(s, s_mean, s_sd);

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), a(loc_ind(i)), log_b(loc_ind(i)), s, reparam_s);
  }

  // ------------- Output z -----------------------
  vector<Type> z(loc_ind.size());
  Type p = 0.1;
  for (int i=0; i<y.size();i++){
    z[i] = a(loc_ind(i))-exp(log_b(loc_ind(i)))/s*(1-pow(-log(1-p), -s));
  }
  ADREPORT(z);

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif


