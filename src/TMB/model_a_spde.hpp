#ifndef model_a_spde_hpp
#define model_a_spde_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a_spde(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Matern)
  */ 
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
 
  // data inputs
  DATA_VECTOR(y); // response vector
  DATA_IVECTOR(n_obs); // number of observations per location
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_IVECTOR(meshidxloc); // indices of the locations in the mesh matrix
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained   to positive: 1 , constrained to be negative: 2, or unconstrained: 3
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov.
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  // parameter list
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER(log_b); // log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER(log_sigma_a); // hyperparameter for the Matern SPDE for a
  PARAMETER(log_kappa_a); // as above

  int n = n_obs.size(); // number of locations
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  
  Type nll = Type(0.0); 
  
  // spde approx
  SparseMatrix<Type> Q = Q_spde(spde, kappa_a);
  Type sigma_marg_a = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_a, 2*nu)); // marginal variance for a
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll = SCALE(GMRF(Q), sigma_a/sigma_marg_a)(mu_a);
  
  // calculate the negative log likelihood
  int start_ind = 0; // index of the first observation of location i in n_obs
  int end_ind = 0; // index of the last observation of location i in n_obs
  if (reparam_s == 0){ // this is the case we are using Gumbel distribution
    for(int i=0;i<n;i++) {
      end_ind += n_obs[i];
      for(int j=start_ind;j<end_ind;j++){
        nll -= gumbel_lpdf<Type>(y[j], a[meshidxloc[i]], log_b);
      }
      start_ind += n_obs[i];
    }
  } else{ // the case where we are using GEV distribution with nonzerio shape parameter
    if (s_sd<9999){ // put a prior on s, or log(s), or log(|s|)
      nll -= dnorm(s, s_mean, s_sd, true);
    }
    if (reparam_s == 1){ // if we have stated that s is constrained to be positive, this implies tha    t we are optimizing log(s)
      s = exp(s);
    } else if (reparam_s == 2){ // if we have stated that s is constrained to be negative, this impl    ies that we are optimizing log(-s)
      s = -exp(s);
    } // if we don't use any reparametrization, then s is unconstrained
    for(int i=0;i<n;i++) {
      end_ind += n_obs[i];
      for(int j=start_ind;j<end_ind;j++){
        nll -= gev_lpdf<Type>(y[j], a[meshidxloc[i]], log_b, s);
      }
      start_ind += n_obs[i];
    }
  }

  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
