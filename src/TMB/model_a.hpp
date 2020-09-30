#ifndef model_a_hpp
#define model_a_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a(objective_function<Type>* obj) {
  /*
   Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: log_tau ~ GP(0, Sigma(sigma, ell))
  */ 
  // data inputs
  DATA_INTEGER(n); // number of observations 
  DATA_VECTOR(y); // response vector: mws
  DATA_SCALAR(y_min); // a random number between 0 and min(y)
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.
  
  PARAMETER_VECTOR(logit_tau); // random effect to be integrated out: transformed location parameters of the GEV model  
  PARAMETER(log_b); // log-transformed scale parameter of the GEV model 
  PARAMETER(log_s); // log-transformed shape parameter of the GEV model
  PARAMETER(log_sigma); // log-transformed amplitude parameter (scalar) of the exponential covariance function
  PARAMETER(log_ell); // log-transformed smoothness parameter (scalar) of the exponential covariance function
  
  Type b = exp(log_b);
  Type s = exp(log_s);
  vector<Type> a = invlogit(logit_tau) * y_min + s/b;
  Type sigma = exp(log_sigma);
  Type ell = exp(log_ell);
  
  
  using namespace density;
  matrix<Type> cov(n,n);
  cov_expo2<Type>(cov, dd, sigma, ell, sp_thres);
  // if (sp_thres == 0){
  //   matrix<Type> tmp = -dd/ell;
  //   matrix<Type> dd_exp = exp(tmp.array());
  //   cov = sigma*dd_exp;
  // } else {
  //   // construct the covariance matrix
  //   int i,j; 
  //   for (i = 0; i < n; i++){
  //     cov(i,i) = sigma;
  //     for (j = 0; j < i; j++){
  //       Type h = dd(i,j);
  //       if (h >= sp_thres) {
  //         cov(i,j) = 0;
  //         cov(j,i) = 0;
  //       } else {
  //         cov(i,j) = sigma*exp(-h/ell);  
  //         cov(j,i) = cov(i,j);
  //       }
  //     }
  //   }
  // }
  
  // Type nll = n*log_b;
  Type nll = Type(0.0);
  // the negloglik of GEV
  for(int i=0;i<n;i++) {
    // Type t = 1 + s * (y[i] - a[i]) / b;
    // nll += pow(t, -1/s) + (s + 1)/s * log(t) ;
    nll -= gev_lpdf<Type>(y[i], a[i], log_b, s);
  }
  
  // the negloglik of location parameter ~ MVNORM
  MVNORM_t<Type> mvn_nll(cov);
  nll += mvn_nll(logit_tau);
  
  return nll;  
   
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
