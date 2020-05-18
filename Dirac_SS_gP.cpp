#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


double rcpp_ratio(const ivec &delta_vec, const int &delta_index, const mat &X, const vec &y, double sigma_sq, const double &g) {
  // initial stuff
  double log_output;
  double out;
  // generate with input
  ivec d_with = delta_vec;
  d_with(delta_index) = 1;
  int k_with = sum(d_with);
  uvec col_id_with(k_with);
  int j = 0;
  for (size_t i = 0; i < delta_vec.n_elem; i++) {
    if (d_with(i) == 1) {
      col_id_with(j) = i;
      j++;
    }
  }
  mat X_with = X.cols(col_id_with);
  // generate without input
  ivec d_wout = delta_vec;
  d_wout(delta_index) = 0;
  int k_wout = sum(d_wout);
  uvec col_id_wout(k_wout);
  int h = 0;
  for (size_t i = 0; i < delta_vec.n_elem; i++) {
    if (d_wout(i) == 1) {
      col_id_wout(h) = i;
      h++;
    }
  }
  mat X_wout = X.cols(col_id_wout);
  // parameters with
  mat inv_A_N_with = (g + 1) / g * X_with.t() * X_with;
  mat A_N_with = inv_sympd(inv_A_N_with);
  vec a_N_with = A_N_with * (X_with.t() * y);
  // parameters without
  if (k_wout == 0) {
    log_output = - as_scalar(a_N_with.t() * inv_A_N_with * a_N_with) / sigma_sq + log(g + 1);
  } else {
    mat inv_A_N_wout = (g + 1) / g * X_wout.t() * X_wout;
    mat A_N_wout = inv_sympd(inv_A_N_wout);
    vec a_N_wout = A_N_wout * (X_wout.t() * y);
    log_output = - as_scalar(a_N_with.t() * inv_A_N_with * a_N_with - a_N_wout.t() * inv_A_N_wout * a_N_wout) / sigma_sq + log(g + 1);
  }
  // output
  out = exp(log_output / 2);
  return out;
}

double r_beta(const double &v, const double &V) {
  Function r_beta("rbeta");
  NumericVector out = r_beta(1, _["shape1"] = v, _["shape2"] = V);
  return out[0];
}


// [[Rcpp::export]]
List rcpp_Dirac_SS_g(const mat &X, const vec &y, const int &n_samples, const double &burn_in,
  const double s_0 = .001, const double S_0 = .001,
  const bool update_psi = 1, const double fixed_psi = 1, const double b_0 = .5, const double B_0 = .5,
  const bool update_omega = 1, const double fixed_omega = .5, const double v_0 = 1, const double V_0 = 1
) {
  // initialization of utilities
  int total_samples = n_samples + burn_in;
  int k = X.n_cols;
  uvec k_perm(k), active_col_id;
  int n = y.n_elem;
  int active_size, j;
  mat X_active(X.n_rows, X.n_cols);
  double pip, s_N, S_N, b_N, B_N, v_N, V_N;
  // output variables
  List out(5);
  imat delta_out(n_samples, k, fill::zeros);
  mat beta_out(n_samples, k, fill::zeros);
  vec sigma_sq_out(n_samples, fill::zeros);
  vec psi_out(n_samples, fill::zeros);
  vec omega_out(n_samples, fill::zeros);
  int saving_iteration = 0;
  // initial guesses
  ivec delta_vec = randi(k, distr_param(0, 1));
  double sigma_sq_draw = var(y);
  double psi_draw;
  mat inv_A_N, A_N;
  vec inv_A_0_v(k);
  vec a_N(k, fill::zeros);
  vec beta_draw(k, fill::zeros);
  vec beta_active_draw(k, fill::zeros);
  if (update_psi) {
    psi_draw = 1 / randg(distr_param(b_0, 1 / B_0));
  } else {
    psi_draw = fixed_psi;
  }
  double omega_draw;
  if (update_omega) {
    omega_draw = r_beta(v_0, V_0);
  } else {
    omega_draw = fixed_omega;
  }
  // start sampler
  for (int sample = 0; sample < total_samples; sample++) {
    // sample active coefficients
    k_perm = randperm(k);
    for (int i = 0; i < k; i++) {
      pip = 1 / (1 + rcpp_ratio(delta_vec, k_perm(i), X, y, sigma_sq_draw, psi_draw) * (1 - fixed_omega) / fixed_omega);
      delta_vec(k_perm(i)) = randu() < pip;
    }
    active_size = sum(delta_vec);
    // sample all other parameters
    if (active_size == 0) {
      s_N = s_0 + (n - 1) / 2;
      S_N = S_0 + dot(y, y) / 2;
      sigma_sq_draw = 1 / randg(distr_param(s_N, 1 / S_N));
    } else {
      // active predictors
      X_active.set_size(X.n_rows, active_size);
      active_col_id.set_size(active_size);
      j = 0;
      for (size_t i = 0; i < k; i++) {
        if (delta_vec(i) == 1) {
          active_col_id(j) = i;
          j++;
        }
      }
      X_active = X.cols(active_col_id);
      // resize beta parameters
      inv_A_0_v.set_size(active_size);
      inv_A_N.set_size(active_size, active_size);
      A_N.set_size(active_size, active_size);
      a_N.set_size(active_size);
      beta_active_draw.set_size(active_size);
      beta_draw.zeros();
      // update beta_draw and beta_active_draw
      inv_A_0_v.fill(1 / psi_draw);
      inv_A_N = X_active.t() * X_active + diagmat(inv_A_0_v);
      A_N = inv_sympd(inv_A_N);
      a_N = A_N * (X_active.t() * y);
      beta_active_draw = mvnrnd(a_N, sigma_sq_draw * A_N);
      for (size_t i = 0; i < active_size; i++) {
        beta_draw(active_col_id(i)) = beta_active_draw(i);
      }
      // update sigma_sq_draw
      s_N = s_0 + (n - 1) / 2;
      S_N = S_0 + (dot(y, y) - as_scalar(a_N.t() * inv_A_N * a_N)) / 2;
      sigma_sq_draw = 1 / randg(distr_param(s_N, 1 / S_N));
    }
    // update psi_draw
    if (update_psi) {
      b_N = b_0 + active_size / 2;
      B_N = B_0 + dot(beta_active_draw, beta_active_draw) / (2 * sigma_sq_draw);
      psi_draw = 1 / randg(distr_param(b_N, 1 / B_N));
    }
    // update omega_draw
    if (update_omega) {
      v_N = v_0 + active_size;
      V_N = V_0 + k - active_size;
      omega_draw = r_beta(v_N, V_N);
    }
    // record draws
    if (sample >= burn_in) {
      delta_out.row(saving_iteration) = delta_vec.t();
      beta_out.row(saving_iteration) = beta_draw.t();
      sigma_sq_out(saving_iteration) = sigma_sq_draw;
      if (update_psi) {psi_out(saving_iteration) = psi_draw;}
      if (update_omega) {omega_out(saving_iteration) = omega_draw;}
      saving_iteration++;
    }
  }
  // prepare and return output
  out(0) = delta_out;
  out(1) = beta_out;
  out(2) = sigma_sq_out;
  if (update_psi) {
    out(3) = psi_out;
  } else {
    out(3) = fixed_psi;
  }
  if (update_omega) {
    out(4) = omega_out;
  } else {
    out(4) = fixed_omega;
  }
  return out;
}
