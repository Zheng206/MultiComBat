data {
  int<lower=1> N;          // Total observations
  int<lower=1> G;          // Number of features
  int<lower=1> I;          // Number of batches
  array[N] real y;         // Observed values
  array[N] int<lower=1> g; // Feature indices
  array[N] int<lower=1> i; // Batch indices

  // Half-t hyperparameters
  real<lower=0> nu_psi;    // df for batch-level
  real<lower=0> nu_tau;    // df for feature-level
  real<lower=0> nu_delta;  // df for observation-level
  real<lower=0> A_psi;     // scale for batch-level
  real<lower=0> A_tau;     // scale for feature-level
  real<lower=0> A_delta;   // scale for observation-level
}

parameters {
  // Means
  vector[I] gamma_i;               // Batch-level means
  matrix[G, I] mu_ig;              // Feature-level means
  
  // Standard deviations (Half-t parameterization)
  vector<lower=0>[I] sigma_psi;         // Batch-level SD
  matrix<lower=0>[G, I] sigma_tau;      // Feature-level SD
  matrix<lower=0>[G, I] sigma_delta;    // Observation-level SD
}

model {
  // Half-t priors on standard deviations
  sigma_psi ~ student_t(nu_psi, 0, A_psi);
  to_vector(sigma_tau) ~ student_t(nu_tau, 0, A_tau);
  to_vector(sigma_delta) ~ student_t(nu_delta, 0, A_delta);

  // Weakly informative prior for batch means
  gamma_i ~ normal(0, 5);

  // Hierarchical structure - simpler without measurement level
  for (gg in 1:G) {
    for (ii in 1:I) {
      mu_ig[gg, ii] ~ normal(gamma_i[ii], sigma_psi[ii]);
    }
  }

  // Likelihood - now directly using mu_ig since there's no measurement level
  for (n in 1:N) {
    y[n] ~ normal(mu_ig[g[n], i[n]], sigma_delta[g[n], i[n]]);
  }
}

generated quantities {
  array[N] real y_rep;
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu_ig[g[n], i[n]], sigma_delta[g[n], i[n]]);
  }
}
