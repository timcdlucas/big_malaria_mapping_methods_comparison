data {
  int<lower=1> n;
  int<lower=1> k;
  int<lower=1> m;
  matrix[n,m] x;
  vector[n] y;
  matrix[k,m] omega;
  real<lower=0> l1;
  real<lower = 0> scale_sd;
}
transformed data {
  
}
parameters {
  vector[k] beta1;
  vector[k] beta2;
  real<lower=0> sigma2;
  vector[m] logbw;
}
transformed parameters {
  vector[n] fhat;
  matrix[n,k] cosfeatures;
  matrix[n,k] sinfeatures;
  matrix[n, m] scale_x;
  real scale;
  vector[m] bw;
  matrix[n,k] features;

  bw = exp(logbw);
  
  for(i in 1:m) {
    scale_x[, i] = x[, i] * bw[i];
  }
   
  features = x * omega';
  
  scale = sqrt(2.0/n);
  for(i in 1:n)
    for(j in 1:k) {
      cosfeatures[i,j] = cos(features[i,j]);
      sinfeatures[i,j] = sin(features[i,j]);
    }
  cosfeatures = cosfeatures * scale;
  sinfeatures = sinfeatures * scale;
  
  fhat = cosfeatures * beta1 + sinfeatures * beta2;
}
model {
  target += normal_lpdf(logbw | 0, scale_sd);
  target += normal_lpdf(beta1 | 0, l1);
  target += normal_lpdf(beta2 | 0, l1);
  target += normal_lpdf(sigma2 | 0, 1);
  y ~ normal(fhat, sigma2);
}
