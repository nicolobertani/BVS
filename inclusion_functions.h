#ifndef _INCLUSION_FNS
#define _INCLUSION_FNS

#include <armadillo>

// independent prior inclusion ratio
double IP_ratio(const ivec &delta_vec, const int &delta_index, const mat &X, const vec &y, double sigma_sq, const double &psi) {
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
  vec inv_A_0_with_v(k_with);
  inv_A_0_with_v.fill(1 / psi);
  mat inv_A_0_with = diagmat(inv_A_0_with_v);
  mat inv_A_N_with = X_with.t() * X_with + inv_A_0_with;
  mat A_N_with = inv_sympd(inv_A_N_with);
  vec a_N_with = A_N_with * (X_with.t() * y);
  // parameters without
  if (k_wout == 0) {
    log_output = - as_scalar(a_N_with.t() * inv_A_N_with * a_N_with) / sigma_sq + log(psi / det(A_N_with));
  } else {
    vec inv_A_0_wout_v(k_wout);
    inv_A_0_wout_v.fill(1 / psi);
    mat inv_A_0_wout = diagmat(inv_A_0_wout_v);
    mat inv_A_N_wout = X_wout.t() * X_wout + inv_A_0_wout;
    mat A_N_wout = inv_sympd(inv_A_N_wout);
    vec a_N_wout = A_N_wout * (X_wout.t() * y);
    log_output = - as_scalar(a_N_with.t() * inv_A_N_with * a_N_with - a_N_wout.t() * inv_A_N_wout * a_N_wout) / sigma_sq + log(psi * det(A_N_wout) / det(A_N_with));
  }
  // output
  out = exp(log_output / 2);
  return out;
}

// g prior inclusion ratio
double gP_ratio(
  const arma::ivec &delta_vec, const int &delta_index,
  const arma::mat &X, const arma::vec &y,
  double sigma_sq, const double &g
) {
  // initial stuff
  double log_output, out;
  // generate with input
  arma::ivec d_with = delta_vec;
  d_with(delta_index) = 1;
  int k_with = arma::sum(d_with);
  arma::uvec col_id_with(k_with);
  int j = 0;
  for (size_t i = 0; i < delta_vec.n_elem; i++) {
    if (d_with(i) == 1) {
      col_id_with(j) = i;
      j++;
    }
  }
  arma::mat X_with = X.cols(col_id_with);
  // generate without input
  arma::ivec d_wout = delta_vec;
  d_wout(delta_index) = 0;
  int k_wout = arma::sum(d_wout);
  arma::uvec col_id_wout(k_wout);
  int h = 0;
  for (size_t i = 0; i < delta_vec.n_elem; i++) {
    if (d_wout(i) == 1) {
      col_id_wout(h) = i;
      h++;
    }
  }
  arma::mat X_wout = X.cols(col_id_wout);
  // parameters with
  arma::mat inv_A_N_with = (g + 1) / g * X_with.t() * X_with;
  arma::mat A_N_with = arma::inv_sympd(inv_A_N_with);
  arma::vec a_N_with = A_N_with * (X_with.t() * y);
  // parameters without
  if (k_wout == 0) {
    log_output =
    - arma::as_scalar(a_N_with.t() * inv_A_N_with * a_N_with) / sigma_sq
    + std::log(g + 1);
  } else {
    arma::mat inv_A_N_wout = (g + 1) / g * X_wout.t() * X_wout;
    arma::mat A_N_wout = inv_sympd(inv_A_N_wout);
    arma::vec a_N_wout = A_N_wout * (X_wout.t() * y);
    log_output =
    - arma::as_scalar(a_N_with.t() * inv_A_N_with * a_N_with - a_N_wout.t() * inv_A_N_wout * a_N_wout) / sigma_sq
    + std::log(g + 1);
  }
  // output
  out = std::exp(log_output / 2);
  return out;
}

#endif
