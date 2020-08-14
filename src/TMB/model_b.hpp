#ifndef model_b_hpp
#define model_b_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_b(objective_function<Type>* obj) {
  /*
   Model layer 1: y ~ GEV(a, b, s)
   Model layer 2: log_b ~ GP(0, Sigma_b(sigma, ell))
   */ 
  // data inputs
  DATA_INTEGER(n); // number of observations 
  DATA_VECTOR(y); // response vector: mws
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.
  
  PARAMETER(a); // location parameters of the GEV model  
  PARAMETER_VECTOR(log_b); // log-transformed scale parameter of the GEV model: the random effect to be integrated out 
  PARAMETER(log_s); // log-transformed shape parameter of the GEV model
  PARAMETER(log_sigma); // log-transformed amplitude parameter (scalar) of the exponential covariance function
  PARAMETER(log_ell); // log-transformed smoothness parameter (scalar) of the exponential covariance function
  
  vector<Type> b = exp(log_b);
  Type s = exp(log_s);
  Type sigma = exp(log_sigma);
  Type ell = exp(log_ell);
  
  
  using namespace density;
 
  matrix<Type> cov(n,n);
  if (sp_thres == 0){
    matrix<Type> tmp = -dd/ell;
    matrix<Type> dd_exp = exp(tmp.array());
    cov = sigma*dd_exp;
  } else {
    // construct the covariance matrix
    int i,j; 
    for (i = 0; i < n; i++){
      cov(i,i) = sigma;
      for (j = 0; j < i; j++){
        Type h = dd(i,j);
        if (h >= sp_thres) {
          cov(i,j) = 0;
          cov(j,i) = 0;
        } else {
          cov(i,j) = sigma*exp(-h/ell);  
          cov(j,i) = cov(i,j);
        }
      }
    }
  }
  
  Type nll = sum(log_b); 
  // the negloglik of GEV
  for(int i=0;i<n;i++) {
    Type t = 1 + s * (y[i] - a) / b[i];
    nll += pow(t, -1/s) + (s + 1)/s * log(t) ;
  }
  
  // the negloglik of scale parameter ~ MVNORM
  MVNORM_t<Type> mvn_nll(cov);
  nll += mvn_nll(log_b);
  
  return nll;  
   
}


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
