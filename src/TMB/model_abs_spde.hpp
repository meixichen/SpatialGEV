#ifndef model_abs_spde_hpp
#define model_abs_spde_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_spde(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: 
  a ~ GP(0, Matern_SPDE)
  logb ~ GP(0, Matern_SPDE)
  g(s) ~ GP(0, Matern_SPDE) where s is a transformation function of s
  */ 
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(n_obs); // number of observations per location
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_MATRIX(design_mat_s); // n x r design matrix for s
  DATA_IVECTOR(meshidxloc); // indices of the locations in the mesh matrix
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained to positive: 1 , constrained to be negative: 2, or unconstrained: 3  
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. 
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_b_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_s_prior); // length 2 vector containing mean and sd of normal prior on beta
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER_VECTOR(s); // If reparam_s = "negative" or "positive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER_VECTOR(beta_s); // r x 1 mean vector coefficients for s
  PARAMETER(log_sigma_a); // hyperparameter for the Matern SPDE for a
  PARAMETER(log_kappa_a); // as above
  PARAMETER(log_sigma_b); // hyperparameter for the Matern SPDE for log(b)
  PARAMETER(log_kappa_b); // as above
  PARAMETER(log_sigma_s); // hyperparameter for the Matern SPDE for s
  PARAMETER(log_kappa_s); // as above

  int n = n_obs.size(); // number of locations
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b); 
  Type kappa_b = exp(log_kappa_b); 
  Type sigma_s = exp(log_sigma_s); 
  Type kappa_s = exp(log_kappa_s); 
 
  Type nll = Type(0.0);
  
  // spde approx
  SparseMatrix<Type> Q_a = Q_spde(spde, kappa_a);
  SparseMatrix<Type> Q_b = Q_spde(spde, kappa_b);
  SparseMatrix<Type> Q_s = Q_spde(spde, kappa_s);
  Type sigma_marg_a = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_a, 2*nu)); // marginal variance for a
  Type sigma_marg_b = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_b, 2*nu)); // marginal variance for log(b)
  Type sigma_marg_s = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_s, 2*nu)); // marginal variance for s
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll = SCALE(GMRF(Q_a), sigma_a/sigma_marg_a)(mu_a);
  nll += SCALE(GMRF(Q_b), sigma_b/sigma_marg_b)(mu_b);
  nll += SCALE(GMRF(Q_s), sigma_s/sigma_marg_s)(mu_s);

  // calculate the negative log likelihood
  int start_ind = 0; // index of the first observation of location i in n_obs
  int end_ind = 0; // index of the last observation of location i in n_obs
  if (reparam_s == 0){ // this is the case we are using Gumbel distribution
    for(int i=0;i<n;i++) {
      end_ind += n_obs[i];
      for(int j=start_ind;j<end_ind;j++){
        nll -= gumbel_lpdf<Type>(y[j], a[meshidxloc[i]], log_b[meshidxloc[i]]);
      }
      start_ind += n_obs[i];
    }
  } else{ // the case where we are using GEV distribution with nonzero shape parameter
    if (reparam_s == 1){ // if we have stated that s is constrained to be positive, this implies that we are optimizing log(s)
      s = exp(s);
    } else if (reparam_s == 2){ // if we have stated that s is constrained to be negative, this implies that we are optimizing log(-s)
      s = -exp(s);
    } // if we don't use any reparametrization, then s is unconstrained
    for(int i=0;i<n;i++) {
      end_ind += n_obs[i];
      for(int j=start_ind;j<end_ind;j++){
        nll -= gev_lpdf<Type>(y[j], a[meshidxloc[i]], log_b[meshidxloc[i]], s[meshidxloc[i]]);
      }
      start_ind += n_obs[i];
    }
  }
  
  // prior
  nll_accumulator_beta<Type>(nll, beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll_accumulator_beta<Type>(nll, beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll_accumulator_beta<Type>(nll, beta_s, beta_prior, beta_s_prior[0], beta_s_prior[1]);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

