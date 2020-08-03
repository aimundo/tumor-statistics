data {
  int<lower=0> N; //total number of observations
  int<lower=0> n_time; //number of repeated measurements
  int<lower=0> n_treatment; //number of treatments
  int<lower=0> q;// wiggliness of the response
  int<lower=0> time_idx[N]; #all timepoints
  int<lower=0> treatment_idx[N];#all observations per group
  vector[N] y; #vector of generated responses
  matrix[n_time, q] X; //polynomial spline matrix
}

parameters {  //parameters to be estimated
  matrix[q, n_treatment] beta;//coefficients??
  real<lower=0> sigma;
}

transformed parameters {
  matrix[n_time, n_treatment] mu = X * beta; //mean response?

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  to_vector(beta) ~ normal(0, 1);
  sigma ~ cauchy(0, 2);
  for (i in 1:N) {
      y[i] ~ normal(mu[time_idx[i], treatment_idx[i]], sigma);
  }
}

generated quantities {
  vector[N] y_rep;
  for (i in 1:N) {
      y_rep[i] = normal_rng(mu[time_idx[i], treatment_idx[i]], sigma);
  }
}

