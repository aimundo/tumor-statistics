data {
  int<lower=0> N;
  int<lower=0> n_time;
  int<lower=0> n_treatment;
  int<lower=0> q;
  int<lower=0> time_idx[N];
  int<lower=0> treatment_idx[N];
  vector[N] y;
  matrix[n_time, q] X;
}

parameters {
  matrix[q, n_treatment] beta;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[n_time, n_treatment] mu = X * beta;

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

