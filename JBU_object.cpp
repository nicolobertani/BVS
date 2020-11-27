#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// r drawing functions
vec r_binom(const int &n, const vec &p) {
  Function r_binom("rbinom");
  NumericVector tmp = r_binom(_["n"] = n, _["size"] = 1, _["prob"] = p);
  vec out(p.n_elem);
  for (uword i = 0; i < p.n_elem; i++) {
    out(i) = tmp[i];
  }
  return out;
}

vec r_gamma(const int &n, const double &shape, const vec &rate) {
  Function r_gamma("rgamma");
  NumericVector tmp = r_gamma(_["n"] = n, _["shape"] = shape, _["rate"] = rate);
  vec out(tmp.n_elem);
  for (uword i = 0; i < tmp.n_elem; i++) {
    out(i) = tmp[i];
  }
  return out;
}


class JBU_model {

  mat m_X_block, m_y_mx;
  umat m_delta;
  sp_mat m_X_alpha;
  // hyperparameters
  double m_psi = pow(10, 4); // beta
  double m_a_1 = 5; // tau
  double m_a_2 = 50; // tau
  double m_epsilon = .005; // gamma
  double m_w = .5; // omega
  double m_v_0;
  mat m_S_0;
  // sampler draws
  vec m_kappa_draw, m_inv_tau_sq_draw, m_beta_draw;
  mat m_inv_Sigma, m_GRR_mx;

  // FILTER
  /* calculate individual approximate posterior probability
  to be used in the loop in the filter */
  double calculate_pip (const uword &i, const uword &j, const double &yTy) {
    double inv_A_N, A_N, a_N, B;
    inv_A_N = dot(m_X_block.col(j), m_X_block.col(j)) + 1 / m_psi;
    A_N = 1 / inv_A_N;
    a_N = A_N * dot(m_X_block.col(j), m_y_mx.col(i));
    B = yTy - pow(a_N, 2) * inv_A_N;
    return log(A_N) - (m_X_block.n_rows - 1) * log(B);
  }
  // SAMPLER
  void set_Wishart_hp (const bool &RATS) {
    if (RATS) {
      m_v_0 = m_y_mx.n_cols * (m_delta.n_rows + 1) - 2; //RATS
    } else {
      m_v_0 = m_y_mx.n_cols + 1; // Jeffrey's
    }
    m_S_0.eye(m_y_mx.n_cols, m_y_mx.n_cols);
    m_S_0 = m_S_0 * m_v_0;
  }
  void set_initial_draws () {
    m_kappa_draw.set_size(m_delta.n_elem);
    m_kappa_draw.ones();
    m_inv_tau_sq_draw.set_size(m_delta.n_elem);
    m_inv_tau_sq_draw.fill(1 / pow(10, 6));
    m_inv_Sigma.set_size(m_y_mx.n_cols, m_y_mx.n_cols);
    m_inv_Sigma.eye();
  }
  void beta_update(const mat &chol_U_Inv_Sigma, const vec &y_vec, const int &blocks) {
    mat UX = chol_U_Inv_Sigma * m_X_alpha;
    vec U_y_vec = chol_U_Inv_Sigma * y_vec;
    for (uword b = 0; b < blocks; b++) {
      uvec selector = ;
      /* code */
    }
  }
  void kappa_update () {
    vec beta_sq, sd_term, w_1, w_2, v;
    beta_sq = pow(m_beta_draw, 2);
    sd_term = .5 * m_inv_tau_sq_draw * beta_sq;
    w_1 = exp(log(1 - m_w) - .5 * log(m_epsilon) - sd_term / m_epsilon);
    w_2 = exp(log(m_w) - sd_term);
    // draw
    /*v = r_binom(m_beta_draw.n_elem, w_1 / (w_1 + w_2));
    m_kappa_draw = m_epsilon * v + (1 - v);*/
    m_kappa_draw = m_beta_draw.n_elem, w_1 / (w_1 + w_2);
  }
  void inv_tau_sq_update () {
  vec beta_sq, a_N_2;
  beta_sq = pow(m_beta_draw, 2);
  a_N_2 = m_a_2 + .5 * beta_sq / m_kappa_draw;
  // draw
  // m_inv_tau_sq_draw = r_gamma(m_inv_tau_sq_draw.n_elem, m_a_1 + .5, a_N_2);
  m_inv_tau_sq_draw = (m_a_1 + .5) / a_N_2;
}
  void GRR_update() { m_GRR_mx = diagmat(m_inv_tau_sq_draw / m_kappa_draw); }
  void w_update() {}
  void inv_Sigma_update () {}




public:
  // FILTER
  // construct and print X and y
  void set_X_block (mat value) { m_X_block = value ;}
  mat get_X_block () { return m_X_block; }
  void set_y_mx (mat value) { m_y_mx = value ;}
  mat get_y_mx () { return m_y_mx; }

  // approximate Dirac filter: populate m_delta
  void filter (const int &K) {
    m_delta.set_size(K, m_y_mx.n_cols);
    vec pip(m_X_block.n_cols, fill::zeros);
    for (uword i = 0; i < m_y_mx.n_cols; i++) {
      double yTy = dot(m_y_mx.col(i), m_y_mx.col(i));
      uvec order_v;
      for (uword j = 0; j < m_X_block.n_cols; j++) {
        pip(j) = calculate_pip(i, j, yTy);
      }
      order_v = sort_index(pip, "descent");
      m_delta.col(i) = sort(order_v.head(K)); // sorting
    }
  }
  umat get_delta_mx () { return m_delta; }

  // create block diagional covariate matrix X_alpha
  void fill_X_alpha () {
    m_X_alpha.set_size(m_X_block.n_rows * m_y_mx.n_cols, m_delta.n_rows * m_y_mx.n_cols);
    for (uword i = 0; i < m_y_mx.n_cols; i++) {
      m_X_alpha.submat(
        i * m_X_block.n_rows, i * m_delta.n_rows,
        (i + 1) * m_X_block.n_rows - 1, (i + 1) * m_delta.n_rows - 1
      ) = m_X_block.cols(m_delta.col(i));
    }
  }
  mat get_X_alpha () { return mat(m_X_alpha); }

  // SAMPLER
  List sample(const bool &RATS, const bool &large, const int &iterations) {
    List out(4);
    // helpers
    mat Kron_hlpr(m_X_block.n_rows, m_X_block.n_rows, fill::eye);
    // hyperparameters Wishart
    set_Wishart_hp(RATS);
    // initial draws
    set_initial_draws();
    GRR_update();
    mat chol_U_Inv_Sigma = kron(m_inv_Sigma, Kron_hlpr);
    vec y_vec = vectorise(m_y_mx);

    for (size_t iter = 0; iter < iterations; iter++) {
      kappa_update();
    }
    out(0) = m_beta_draw;
    out(1) = m_kappa_draw;
    out(2) = m_inv_tau_sq_draw;
    out(3) = m_inv_Sigma;
    return out;
  }

};

/*
END OF JBU MODEL CLASS
START R FUNCTIONS
*/

// [[Rcpp::export]]
umat JBU_filter (const mat &X_block, const mat &y_mx, const int &K) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  return mod.get_delta_mx() + 1;
}

// [[Rcpp::export]]
mat JBU_X_alpha(const mat &X_block, const mat &y_mx, const int &K) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  mod.fill_X_alpha();
  return mod.get_X_alpha();
}

// [[Rcpp::export]]
List test(const mat &X_block, const mat &y_mx, const int &K, const bool &RATS, const bool &large) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  mod.fill_X_alpha();
  List out = mod.sample(RATS, large, 1);
  return out;
}
