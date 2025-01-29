functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1) {
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    }
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }

  matrix construct_spline_bases(int N, real[] X, int n_knot, vector knots, int degree);
  matrix construct_spline_bases(int N, real[] X, int n_knot, vector knots, int degree) {

    int n_basis = n_knot + degree - 1; // total number of B-splines

    vector[degree + n_knot] ext_knots_temp;
    vector[2*degree + n_knot] ext_knots; // set of extended knots
    ext_knots_temp = append_row(rep_vector(knots[1], degree), knots);
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[n_knot], degree));

    matrix[n_basis, N] B; // matrix of B-splines
    for (ind in 1:n_basis)
      B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, degree + 1));

    B[n_basis, N] = 1; // what is this line for??
    return B;
  }
}

data {
  // input data
  int<lower=0> N;                 // number of observations
  int<lower=0> K;                 // number of cell types
  matrix<lower=0, upper=1>[N, K] X;           // reference panel
  real<lower=0> z[N];              // CpG density of reference sites
  
  // spline parameters
  int n_knot_z; // num of knots
  vector[n_knot_z] knots_z; // the sequence of knots
  int degree_z; // the degree of spline (is equal to order - 1)
  
  // cell type proportion
  simplex[K] pi;
  
  // regression coefficients
  int<lower=0> L;
  vector[L] w_mu;
  vector[L] w_sigma;

}

transformed data {

  // Spline basis for covariate z (i.e., number of CpGs)
  int n_basis_z = n_knot_z + degree_z - 1;
  // Spline basis for z
  matrix[n_basis_z, N] B_z;
  B_z = construct_spline_bases(N, z, n_knot_z, knots_z, degree_z);

  // // Number of predictors in the design matrix for mu and sigma
  // int L_mu = n_basis_z + 1;
  // int L_sigma = n_basis_z + 1;
}

generated quantities {

  if (L != n_basis_z + 1) {
    reject("L must be n_basis_z + 1"); // This will throw an error if condition fails
  }
  
  vector[N] xtilde = X * pi;

  // Design matrix for location
  matrix[N, L] D_mu;
  for (k in 1:(L-1)) {
    D_mu[:, k] = xtilde .* to_vector(B_z[k, :]);
  }
  D_mu[:, L] = rep_vector(1, N); // the intercept

  // Design matrix for dispersion
  matrix[N, L] D_sigma;
  for (k in 1:(L-1)) {
    D_sigma[:, k] = xtilde .* to_vector(B_z[k, :]);
  }
  D_sigma[:, L] = rep_vector(1, N); // the intercept

  // GLM Likelihood
  vector[N] mu = exp(D_mu * w_mu);
  vector[N] sigma = exp(D_sigma * w_sigma);
  int y_sim[N];

  for (n in 1:N) {
    y_sim[n] = neg_binomial_2_rng(mu[n], sigma[n]);
  }
}
