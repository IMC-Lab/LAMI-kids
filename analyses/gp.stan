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
  array[N] int<lower=0> y;  // output data (counts)
}

transformed data {
  // pre-compute squared distance matrix
  array[G, G] vector[D] dist_sq;
  for (i in 1:G)
    for (j in 1:G)
      for (d in 1:D)
	dist_sq[i, j, d] = square(grid[i, d] - grid[j, d]);
}

parameters {
  real a;                        // global intercept mean
  array[C] vector[G] eta;              // population-level GP mean variates (by condition)
  array[C, P] vector[G] eta_tilde;     // participant-level GP mean variates

  // Hyperparameters
  vector<lower=0>[D] rho;        // population-level length-scale
  real<lower=0> alpha;           // population-level marginal standard deviation
  vector<lower=0>[D] rho_tilde;  // participant-level length-scale
  real<lower=0> alpha_tilde;     // participant-level marginal standard deviation
}

transformed parameters {
  array[C] vector[G] f;                // group-level GPs (by condition)
  array[C, P] vector[G] f_tilde;       // participant-level GPs (by participantXcondition)
  
  {
    // pre-compute GP covariance matrices
    matrix[G, G] Lf = cov_exp_quad_ARD(dist_sq, alpha, rho);
    matrix[G, G] Lf_tilde = cov_exp_quad_ARD(dist_sq, alpha_tilde, rho_tilde);

    for (i in 1:C) {
      f[i] = Lf * eta[i];
      for (j in 1:P) {
        f_tilde[i, j] = Lf_tilde * eta_tilde[i, j];
      }
    }
  }
}

model {
  // intercept terms
  a ~ normal(0, 5);
  
  // GP terms
  for (i in 1:C) {
    eta[i] ~ std_normal();
    for (j in 1:P) {
      eta_tilde[i, j] ~ std_normal();
    }
  }

  // Hyperparameters
  rho ~ inv_gamma(5, 500);
  alpha ~ std_normal();
  rho_tilde ~ inv_gamma(5, 500);
  alpha_tilde ~ std_normal();

  {
    vector[N] lambda;
    for (i in 1:N) {
      lambda[i] = a + f[c[i], g[i]] + f_tilde[c[i], p[i], g[i]];
    }

    y ~ poisson_log(lambda);
  }
}
