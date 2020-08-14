#ifndef model_ab_hpp
#define model_ab_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_ab(objective_function<Type>* obj) {
  /*
   Model layer 1: y ~ GEV(a, b, s)
   Model layer 2: 
   logit_tau ~ GP(0, Sigma_a(sigma_a, ell_a))
   logb ~ GP(0, Sigma_b(sigma_b, ell_b))
   */ 
  // data inputs
  DATA_INTEGER(n); // number of observations 
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_SCALAR(y_min); // a random number between 0 and min(y)
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.

  PARAMETER_VECTOR(logit_tau); // random effect to be integrated out: transformed location parameters of the GEV model  
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER(log_s); // log-transformed shape parameters of the GEV model
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_sigma_b); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_b
  PARAMETER(log_ell_b); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_b
  
  vector<Type> b = exp(log_b);
  Type s = exp(log_s);
  vector<Type> a = invlogit(logit_tau) * y_min + s/b;
  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);
  Type sigma_b = exp(log_sigma_b);
  Type ell_b = exp(log_ell_b);

  using namespace density;
  
  // construct the covariance matrices
  matrix<Type> cova(n,n);
  matrix<Type> covb(n,n);
  matrix<Type> tmp(n,n);
  matrix<Type> dd_exp(n,n);
  if (sp_thres == 0){
    tmp = -dd/ell_a;
    dd_exp = exp(tmp.array());
    cova = sigma_a*dd_exp;
    
    tmp = -dd/ell_b;
    dd_exp = exp(tmp.array());
    covb = sigma_b*dd_exp;
  } else {
    int i,j; 
    for (i = 0; i < n; i++){
      cova(i,i) = sigma_a;
      
      covb(i,i) = sigma_b;
      for (j = 0; j < i; j++){
        Type h = dd(i,j);
        if (h >= sp_thres) {
          cova(i,j) = 0;
          cova(j,i) = 0;
          
          covb(i,j) = 0;
          covb(j,i) = 0;
        } else {
          cova(i,j) = sigma_a*exp(-h/ell_a);  
          cova(j,i) = cova(i,j);
          
          covb(i,j) = sigma_b*exp(-h/ell_b);  
          covb(j,i) = covb(i,j);
        }
      }
    }
  }
  
  Type nll = sum(log_b); 
  // the negloglik of y ~ GEV
  for(int i=0;i<n;i++) {
    Type t = 1 + s * (y[i] - a[i]) / b[i];
    nll += pow(t, -1/s) + (s + 1)/s * log(t) ;
  }
  
  // the negloglik of a ~ GP
  nll += MVNORM(cova)(logit_tau);
  
  // the negloglik of b ~ GP
  nll += MVNORM(covb)(log_b);
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif