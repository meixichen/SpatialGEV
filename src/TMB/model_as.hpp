#ifndef model_as_hpp
#define model_as_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_as(objective_function<Type>* obj) {
    // data inputs
    /*
      Model layer 1: y ~ GEV(a, b, s)
    Model layer 2: 
      logit_tau ~ GP(0, Sigma_b(sigma_b, ell_b))
    logs ~ GP(0, Sigma_s(sigma_s, ell_s)) 
    */ 
    DATA_INTEGER(n); // number of observations 
    DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
    DATA_SCALAR(y_min); // a random number between 0 and min(y)
    DATA_MATRIX(dd); // distance matrix
    DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.
    
    PARAMETER_VECTOR(logit_tau); // random effect to be integrated out: transformed location parameters of the GEV model  
    PARAMETER(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
    PARAMETER_VECTOR(log_s); // random effect to be integrated out: log-transformed shape parameters of the GEV model
    PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
    PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a
    PARAMETER(log_sigma_s); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_s
    PARAMETER(log_ell_s); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_s
    
    Type b = exp(log_b);
    vector<Type> s = exp(log_s);
    vector<Type> a = invlogit(logit_tau) * y_min + s/b;
    Type sigma_a = exp(log_sigma_a);
    Type ell_a = exp(log_ell_a);
    Type sigma_s = exp(log_sigma_s);
    Type ell_s = exp(log_ell_s);  
    
    using namespace density;
    
    // construct the covariance matrices
  matrix<Type> cova(n,n);
  matrix<Type> covs(n,n);
  matrix<Type> tmp(n,n);
  matrix<Type> dd_exp(n,n);
  if (sp_thres == 0){
    tmp = -dd/ell_a;
    dd_exp = exp(tmp.array());
    cova = sigma_a*dd_exp;

    tmp = -dd/ell_s;
    dd_exp = exp(tmp.array());
    covs = sigma_s*dd_exp;
    
  } else {
    int i,j; 
    for (i = 0; i < n; i++){
      cova(i,i) = sigma_a;
      
      covs(i,i) = sigma_s;
      for (j = 0; j < i; j++){
        Type h = dd(i,j);
        if (h >= sp_thres) {
          cova(i,j) = 0;
          cova(j,i) = 0;

          covs(i,j) = 0;
          covs(j,i) = 0;
        } else {
          cova(i,j) = sigma_a*exp(-h/ell_a);  
          cova(j,i) = cova(i,j);
          
          covs(i,j) = sigma_s*exp(-h/ell_s);  
          covs(j,i) = covs(i,j);
        }
      }
    }
  }
    
    Type nll = n*log_b; 
    // the negloglik of y ~ GEV
    for(int i=0;i<n;i++) {
      Type t = 1 + s[i] * (y[i] - a[i]) / b;
      nll += pow(t, -1/s[i]) + (s[i] + 1)/s[i] * log(t) ;
    }
    
    // the negloglik of logit_tau ~ GP
    nll += MVNORM(cova)(logit_tau);
    
    // the negloglik of s ~ GP
    nll += MVNORM(covs)(log_s);
    
    return nll;  
  }

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
