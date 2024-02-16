#ifndef model_a_matern_hpp
#define model_a_matern_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// TMB specification of GEV-GP models with a chosen covariance kernel.
///
/// The model is defined as follows:
///
/// y ~ GEV(a, b, s),
/// a ~ GP(log_sigma_a, log_kappa_a)
/// where the GP is parameterized using the matern covariance kernel.
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
/// @param[in] nu Presepecified smoothness parameter for the Matérn covariance kernel applicable to
/// all random effects.
/// @param[in] a_pc_prior Integer specifying the type of prior to use on the Matérn
/// GP on a.
/// 1 for using PC prior on a, 0 for using Lebesgue prior.
/// `log_sigma_a` and `log_kappa_a`.  1 for using PC prior on a, 0 for using Lebesgue prior.
/// @param[in] range_a_prior PC prior on the range parameter for the Matérn GP on
/// a. Vector of length 2 `(rho_0, p_rho)` s.t. `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_a_prior PC prior on the variance parameter for the Matérn GP on
/// a. Vector of length 2 `(sig_0, p_sig)` s.t. `Pr(sig > sig_0) = p_sig`.
/// @param[in] s_mean Scalar for Normal prior mean on s.
/// @param[in] s_sd Scalar for Normal prior sd on s.
///
/// --------- Parameters to estimate ------------
/// @param[in] a GEV location parameter.
/// Vector of length `n_loc`.
/// @param[in] log_b GEV scale parameter on the log scale.
/// Vector of length 1.
/// @param[in] s GEV shape parameter on the scale specified by `reparam_s`.
/// Vector of length 1.
/// @param[in] beta_a GP mean covariate coefficient vector of length `n_covariate`
/// for a.
/// @param[in] log_sigma_a GP covariance kernel variance hyperparameter
/// for a.
/// @param[in] log_kappa_a GP covariance kernel range hyperparameter
/// for a.
template<class Type>
Type model_a_matern(objective_function<Type>* obj){
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
  DATA_SCALAR(nu);

  // Inputs for a
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);
  DATA_INTEGER(a_pc_prior);
  DATA_VECTOR(range_a_prior);
  DATA_VECTOR(sigma_a_prior);
  DATA_SCALAR(s_mean);
  DATA_SCALAR(s_sd);

  // ------------ Parameters ----------------------

  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);

  PARAMETER_VECTOR(beta_a);
  PARAMETER(log_sigma_a);
  PARAMETER(log_kappa_a);

  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  // ---------- Likelihood contribution from a ------------------
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += nlpdf_gp_matern<Type>(mu_a, dist_mat,
				   exp(log_sigma_a),
				   exp(log_kappa_a),
                                   nu, sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior(0),
                                beta_a_prior(1));
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_a,
					   log_sigma_a,
					   a_pc_prior,
                                           nu, range_a_prior,
					   sigma_a_prior);
  // FIXME: rename this to not depend on `s`
  nll += nlpdf_s_prior<Type>(s(0), s_mean, s_sd);

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), a(loc_ind(i)), log_b(0), s(0), reparam_s);
  }

  // ------------- Output z -----------------------
  vector<Type> z(loc_ind.size());
  Type p = 0.1;
  for (int i=0; i<y.size();i++){
    z[i] = a(loc_ind(i))-exp(log_b(0))/s(0)*(1-pow(-log(1-p), -s(0)));
  }
  ADREPORT(z);

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif


