data {
  int<lower=1> K;  
  int<lower=1> S;  
  int<lower=1> N;  
  int<lower=1> T;  
  array[N, T] int<lower=1, upper=S> y; 
}

parameters {
  
  simplex[K] theta; 
  
  
  array[K, S] simplex[S] A;
}

model {
  
  for (n in 1:N) {
   
    vector[K] lps;
    
    for (k in 1:K) {
      
      lps[k] = log(theta[k]);
      
      
      if (T > 1) {
        for (t in 2:T) {
          lps[k] += categorical_lpmf(y[n, t] | A[k, y[n, t-1]]);
        }
      }
    }
    
    
    target += log_sum_exp(lps);
  }
}