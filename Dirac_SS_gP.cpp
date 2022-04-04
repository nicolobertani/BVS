/*
High performance Bayesian Variable Selection for R using C++ via Rcpp and RcppArmadillo
Copyright (C) 2020  Nicolò Bertani

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "sampling_functions.h"
#include "inclusion_functions.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List Dirac_SS_g(const mat &X, const vec &y, const int &n_samples, const double &burn_in,
  const double s_0 = .001, const double S_0 = .001,
  const bool update_g = 0, const double fixed_g = 1, const double b_0 = .5, const double B_0 = .5,
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
  vec g_out(n_samples, fill::zeros);
  vec omega_out(n_samples, fill::zeros);
  int saving_iteration = 0;
  // initial guesses
  ivec delta_vec = randi(k, distr_param(0, 1));
  double sigma_sq_draw = var(y);
  double g_draw;
  mat inv_A_N, A_N;
  vec a_N(k, fill::zeros);
  vec beta_draw(k, fill::zeros);
  vec beta_active_draw(k, fill::zeros);
  if (update_g) {
    g_draw = 1 / randg(distr_param(b_0, 1 / B_0));
  } else {
    g_draw = fixed_g;
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
      pip = 1 / (1 + gP_ratio(delta_vec, k_perm(i), X, y, sigma_sq_draw, g_draw) * (1 - fixed_omega) / fixed_omega);
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
      inv_A_N.set_size(active_size, active_size);
      A_N.set_size(active_size, active_size);
      a_N.set_size(active_size);
      beta_active_draw.set_size(active_size);
      beta_draw.zeros();
      // update beta_draw and beta_active_draw
      inv_A_N = (g_draw + 1) / g_draw * X_active.t() * X_active;
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
    if (update_g) {
      b_N = b_0 + active_size / 2;
      B_N = B_0 + dot(beta_active_draw, beta_active_draw) / (2 * sigma_sq_draw);
      g_draw = 1 / randg(distr_param(b_N, 1 / B_N));
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
      if (update_g) {g_out(saving_iteration) = g_draw;}
      if (update_omega) {omega_out(saving_iteration) = omega_draw;}
      saving_iteration++;
    }
  }
  // prepare and return output
  out(0) = delta_out;
  out(1) = beta_out;
  out(2) = sigma_sq_out;
  if (update_g) {
    out(3) = g_out;
  } else {
    out(3) = fixed_g;
  }
  if (update_omega) {
    out(4) = omega_out;
  } else {
    out(4) = fixed_omega;
  }
  return out;
}
