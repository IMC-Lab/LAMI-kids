functions {
  /*
   * compute the covariance matrix for a latent Gaussian process
   * args:
   *   dist_sq: the squared distances between GP predictors
   *   alpha: the gp marginal SD parameter
   *   rho: the length-scales for each dimension of grid
   * returns:
   *   a covariance matrix between locations on the grid
   */
  matrix cov_exp_quad_ARD(array[,] vector dist_sq, real alpha, vector rho) {
    int G = size(dist_sq);
    int D = rows(rho);
    real delta = 1e-8;
    matrix[G, G] K;
    real alpha_sq = square(alpha);
    vector[D] rho_sq;
    for (d in 1:D) rho_sq[d] = square(rho[d]);

    for (i in 1:(G-1)) {
      K[i, i] = alpha_sq + delta;
      for (j in (i + 1):G) {
        K[i, j] = alpha_sq * exp(-0.5 * sum(dist_sq[i, j] ./ rho_sq));
        K[j, i] = K[i, j];
      }
    }
    K[G, G] = alpha_sq + delta;
    return cholesky_decompose(K);
  }
}

data {
  int<lower=1> N;     // total number of data points
  int<lower=1> G;     // number of data points per trial (grid size)
  int<lower=1> P;     // number of participants
  int<lower=1> D;     // number of GP input dimensions
  int<lower=1> C;     // number of conditions

  array[N] int<lower=1, upper=G> g;  // grid index for each trial
  array[N] int<lower=1, upper=P> p;  // participant index for each trial
  array[N] int<lower=1, upper=C> c;  // condition index for each trial

  array[G] vector[D] grid;  // input data (same across trials)
}

transformed data {
  // pre-compute squared distance matrix
  array[G, G] vector[D] dist_sq;
  for (i in 1:G)
    for (j in 1:G)
      for (d in 1:D)
	dist_sq[i, j, d] = square(grid[i, d] - grid[j, d]);
}

model {}

generated quantities {
  //vector[N] prior_lambda;
  array[C] vector[G] prior_f;
  array[C, P] vector[G] prior_f_tilde;
  
  // sample from priors
  real prior_a = normal_rng(0, 5);
  vector<lower=0>[D] prior_rho = [inv_gamma_rng(5, 500),
				  inv_gamma_rng(5, 500)]';
  real prior_alpha = normal_rng(0, 1);
  vector<lower=0>[D] prior_rho_tilde = [inv_gamma_rng(5, 500),
					inv_gamma_rng(5, 500)]';
  real prior_alpha_tilde = normal_rng(0, 1);
  {
    // pre-compute GP covariance matrices
    matrix[G, G] prior_Lf = cov_exp_quad_ARD(dist_sq, prior_alpha, prior_rho);
    matrix[G, G] prior_Lf_tilde = cov_exp_quad_ARD(dist_sq, prior_alpha_tilde,
						   prior_rho_tilde);
    array[C] vector[G] prior_eta;
    array[C, P] vector[G] prior_eta_tilde;
    
    for (i in 1:C) {
      for (j in 1:G) prior_eta[i, j] = normal_rng(0, 1);
      prior_f[i] = prior_Lf * prior_eta[i];
      for (j in 1:P) {
	for (k in 1:G) prior_eta_tilde[i, j, k] = normal_rng(0, 1);
        prior_f_tilde[i, j] = prior_Lf_tilde * prior_eta_tilde[i, j];
      }
    }
  }

  // trial-level predictions
  //for (i in 1:N) {
  //  prior_lambda[i] = prior_a + prior_f[c[i], g[i]] + prior_f_tilde[c[i], p[i], g[i]];
  //}
}
