#ifndef model_a_hpp
#define model_a_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Sigma_a(sigma_a, ell_a))
  */ 
  // data inputs
  DATA_INTEGER(n); // number of observations 
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.
  DATA_STRING(reparam_s); // a flag indicating whether the shape parameter is "zero", "unconstrained", constrained to be "negative", or constrained to be "positve"
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER(log_b); // log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a

  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);
  
  using namespace density;
  
  // construct the covariance matrix
  matrix<Type> cov(n,n);
  cov_expo<Type>(cov, dd, sigma_a, ell_a, sp_thres);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  if (reparam_s == "zero"){ // this is the case we are using Gumbel distribution
    for(int i=0;i<n;i++) {
      nll -= gumbel_lpdf<Type>(y[i], a[i], log_b);
    }
  } else{ // the case where we are using GEV distribution with nonzerio shape parameter
    
    if (reparam_s == "positive"){ // if we have stated that s is constrained to be positive, this implies that we are optimizing log(s)
      s = exp(s);
    } else if (reparam_s == "negative"){ // if we have stated that s is constrained to be negative, this implies that we are optimizing log(-s)
      s = -exp(s);
    } // if we don't use any reparametrization, then s is unconstrained
    
    for(int i=0;i<n;i++) {
      nll -= gev_lpdf<Type>(y[i], a[i], log_b, s);
    }
    
  } 
  
  nll += MVNORM(cova)(a);
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
