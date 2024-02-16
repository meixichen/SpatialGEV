#ifndef model_abs_exp_hpp
#define model_abs_exp_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// TMB specification of GEV-GP models with a chosen covariance kernel.
///
/// The model is defined as follows:
///
/// y ~ GEV(a, b, s),
/// a ~ GP(log_sigma_a, log_ell_a)
/// log_b ~ GP(log_sigma_b, log_ell_b)
/// s ~ GP(log_sigma_s, log_ell_s)
/// where the GP is parameterized using the exp covariance kernel.
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
/// @param[in] dist_mat `n_loc x n_loc` distance matrix typically constructed via
/// `stats::dist(coordinates)`.
/// @param[in] sp_thres Scalar number used to make the covariance matrix sparse by thresholding.
/// If sp_thres=-1, no thresholding is made.
/// @param[in] design_mat_a Design matrix of size `n_loc x n_covariate` for parameter
/// a.
/// @param[in] beta_a_prior Vector of length 2 containing the mean and sd of the normal
/// prior on `beta_a`.
/// @param[in] design_mat_b Design matrix of size `n_loc x n_covariate` for parameter
/// log_b.
/// @param[in] beta_b_prior Vector of length 2 containing the mean and sd of the normal
/// prior on `beta_b`.
/// @param[in] design_mat_s Design matrix of size `n_loc x n_covariate` for parameter
/// s.
/// @param[in] beta_s_prior Vector of length 2 containing the mean and sd of the normal
/// prior on `beta_s`.
///
/// --------- Parameters to estimate ------------
/// @param[in] a GEV location parameter.
/// Vector of length `n_loc`.
/// @param[in] log_b GEV scale parameter on the log scale.
/// Vector of length `n_loc`.
/// @param[in] s GEV shape parameter on the scale specified by `reparam_s`.
/// Vector of length `n_loc`.
/// @param[in] beta_a GP mean covariate coefficient vector of length `n_covariate`
/// for a.
/// @param[in] log_sigma_a GP covariance kernel variance hyperparameter
/// for a.
/// @param[in] log_ell_a GP covariance kernel range hyperparameter
/// for a.
/// @param[in] beta_b GP mean covariate coefficient vector of length `n_covariate`
/// for log_b.
/// @param[in] log_sigma_b GP covariance kernel variance hyperparameter
/// for log_b.
/// @param[in] log_ell_b GP covariance kernel range hyperparameter
/// for log_b.
/// @param[in] beta_s GP mean covariate coefficient vector of length `n_covariate`
/// for s.
/// @param[in] log_sigma_s GP covariance kernel variance hyperparameter
/// for s.
/// @param[in] log_ell_s GP covariance kernel range hyperparameter
/// for s.
template<class Type>
Type model_abs_exp(objective_function<Type>* obj){
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;

  // ------ Data inputs ------------
  DATA_VECTOR(y);
  DATA_IVECTOR(loc_ind);
  DATA_INTEGER(reparam_s);
  DATA_INTEGER(beta_prior);
  DATA_MATRIX(dist_mat);
  DATA_SCALAR(sp_thres);

  // Inputs for a
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);
  // Inputs for log_b
  DATA_MATRIX(design_mat_b);
  DATA_VECTOR(beta_b_prior);
  // Inputs for s
  DATA_MATRIX(design_mat_s);
  DATA_VECTOR(beta_s_prior);

  // ------------ Parameters ----------------------

  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);

  PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(beta_b);
  PARAMETER_VECTOR(beta_s);
  PARAMETER(log_sigma_a);
  PARAMETER(log_ell_a);
  PARAMETER(log_sigma_b);
  PARAMETER(log_ell_b);
  PARAMETER(log_sigma_s);
  PARAMETER(log_ell_s);

  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  // ---------- Likelihood contribution from a ------------------
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += nlpdf_gp_exp<Type>(mu_a, dist_mat,
				   exp(log_sigma_a),
				   exp(log_ell_a),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior(0),
                                beta_a_prior(1));
  // ---------- Likelihood contribution from log_b ------------------
  // GP latent layer
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += nlpdf_gp_exp<Type>(mu_b, dist_mat,
				   exp(log_sigma_b),
				   exp(log_ell_b),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior(0),
                                beta_b_prior(1));
  // ---------- Likelihood contribution from s ------------------
  // GP latent layer
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll += nlpdf_gp_exp<Type>(mu_s, dist_mat,
				   exp(log_sigma_s),
				   exp(log_ell_s),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior(0),
                                beta_s_prior(1));

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), a(loc_ind(i)), log_b(loc_ind(i)), s(loc_ind(i)), reparam_s);
  }

  // ------------- Output z -----------------------
  vector<Type> z(loc_ind.size());
  Type p = 0.1;
  for (int i=0; i<y.size();i++){
    z[i] = a(loc_ind(i))-exp(log_b(loc_ind(i)))/s(loc_ind(i))*(1-pow(-log(1-p), -s(loc_ind(i))));
  }
  ADREPORT(z);

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif


