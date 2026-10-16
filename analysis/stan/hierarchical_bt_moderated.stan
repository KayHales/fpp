// Hierarchical Davidson-Bradley-Terry with two additions over hierarchical_bt.stan:
//
//  1. Rater discrimination is predicted by rater covariates, so the model can be
//     asked directly whether some groups separate faces more sharply than others.
//     This is the construct the trans-expertise hypothesis is actually about.
//
//  2. Face scores carry a per-group deviation drawn from N(0, tau). tau is then
//     a posterior answer to "how much do rater groups differ in how they rank
//     faces", which is an upper bound rather than a failure to reject.

data {
  int<lower=1> N;
  int<lower=1> F;
  int<lower=1> R;
  int<lower=1> G;                               // rater groups for the invariance test
  int<lower=0> K;                               // rater covariates predicting discrimination
  array[N] int<lower=1, upper=F> left_face;
  array[N] int<lower=1, upper=F> right_face;
  array[N] int<lower=1, upper=R> rater;
  array[N] int<lower=1, upper=3> y;
  array[R] int<lower=1, upper=G> rater_group;
  matrix[R, K] X;                               // centred covariates
}

parameters {
  vector[F] beta_raw;
  real<lower=0> sigma_beta;
  real theta;
  vector[R] log_gamma_raw;
  real<lower=0> sigma_gamma;
  vector[K] alpha;                              // covariate effects on log discrimination
  real log_nu;
  vector[R] delta_raw;
  real<lower=0> sigma_delta;
  matrix[F, G] u_raw;                           // group deviation per face
  real<lower=0> tau;                            // between-group SD of face scores
}

transformed parameters {
  vector[F] beta = beta_raw * sigma_beta;
  vector[R] log_gamma = (K > 0 ? X * alpha : rep_vector(0.0, R)) + log_gamma_raw * sigma_gamma;
  vector[R] delta = delta_raw * sigma_delta;
  matrix[F, G] u = u_raw * tau;
}

model {
  beta_raw ~ std_normal();
  log_gamma_raw ~ std_normal();
  delta_raw ~ std_normal();
  to_vector(u_raw) ~ std_normal();
  sigma_beta ~ exponential(1);
  sigma_gamma ~ exponential(3);
  sigma_delta ~ exponential(1);
  tau ~ exponential(4);                         // shrinks hard toward invariance
  alpha ~ normal(0, 0.5);
  theta ~ normal(0, 0.5);
  log_nu ~ normal(0, 1);

  {
    array[N] int gidx;
    for (n in 1:N) gidx[n] = rater_group[rater[n]];
    vector[N] g = exp(log_gamma[rater]);
    vector[N] bl; vector[N] br;
    for (n in 1:N) {
      bl[n] = beta[left_face[n]]  + u[left_face[n],  gidx[n]];
      br[n] = beta[right_face[n]] + u[right_face[n], gidx[n]];
    }
    vector[N] a = g .* (bl + theta);
    vector[N] b = g .* br;
    vector[N] t = log_nu + delta[rater] + 0.5 * (a + b);
    for (n in 1:N) y[n] ~ categorical_logit([a[n], b[n], t[n]]');
  }
}
