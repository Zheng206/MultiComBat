data {
  int<lower=1> N;                       // Number of subjects
  int<lower=1> G;                       // Number of features
  int<lower=1> M;                       // Number of measurements
  int<lower=1> I;                       // Number of batches
  array[N] int<lower=1, upper=I> i_id;  // Batch assignment
  array[N, G] vector[M] y;              // Observed data
}

transformed data {
  real<lower=0> nu_alpha = 2.0;         // Gamma shape for nu_sigma
  real<lower=0> nu_beta = 0.1;          // Gamma rate for nu_sigma
}

parameters {
  // Batch-level intercept
  vector[I] gamma_base;                 
  
  // Feature/batch-specific means
  array[I, G] vector[M] mu_ig;
  
  // Standard deviations (Half-t prior)
  array[I, G] vector<lower=1e-3>[M] sigma_ig;
  
  // Correlation matrix (LKJ prior)
  array[I, G] cholesky_factor_corr[M] L_R_ig;
  
  // Hierarchical parameters for sigma_ig
  real<lower=0> nu_sigma;               // Degrees of freedom for Half-t
  vector<lower=0>[M] scale_sigma;       // Scale for Half-t
}

transformed parameters {
  array[I, G] cholesky_factor_cov[M] L_Sigma_ig;
  array[I, G] cov_matrix[M] Sigma_ig;
  
  for (i in 1:I) {
    for (g in 1:G) {
      // Reconstruct covariance matrix
      L_Sigma_ig[i, g] = diag_pre_multiply(sigma_ig[i, g], L_R_ig[i, g]);
      Sigma_ig[i, g] = multiply_lower_tri_self_transpose(L_Sigma_ig[i, g]);
    }
  }
}

model {
  // Hyperpriors
  nu_sigma ~ gamma(nu_alpha, nu_beta);
  scale_sigma ~ cauchy(0, 2.5);         // Weakly informative
  
  // Priors
  gamma_base ~ std_normal();
  
  for (i in 1:I) {
    for (g in 1:G) {
      mu_ig[i, g] ~ normal(rep_vector(gamma_base[i], M), 2);
      
      // Half-t prior for standard deviations
      sigma_ig[i, g] ~ student_t(nu_sigma, 0, scale_sigma);
      
      // LKJ prior for correlation matrix (eta=1 → uniform)
      L_R_ig[i, g] ~ lkj_corr_cholesky(1);
    }
  }
  
  // Likelihood
  for (n in 1:N) {
    int i = i_id[n];
    for (g in 1:G) {
      y[n, g] ~ multi_normal_cholesky(mu_ig[i, g], L_Sigma_ig[i, g]);
    }
  }
}

generated quantities {
  array[N, G] vector[M] y_rep;
  array[I, G] cov_matrix[M] Sigma_ig_out;
  
  for (i in 1:I) {
    for (g in 1:G) {
      Sigma_ig_out[i, g] = Sigma_ig[i, g];
    }
  }
  
  for (n in 1:N) {
    int i = i_id[n];
    for (g in 1:G) {
      y_rep[n, g] = multi_normal_cholesky_rng(mu_ig[i, g], L_Sigma_ig[i, g]);
    }
  }
}


