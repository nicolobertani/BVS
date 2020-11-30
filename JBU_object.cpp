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
  vec out(rate.n_elem);
  for (uword i = 0; i < rate.n_elem; i++) {
    out(i) = tmp[i];
  }
  return out;
}

double r_beta(const double &v, const double &V) {
  Function r_beta("rbeta");
  NumericVector out = r_beta(1, _["shape1"] = v, _["shape2"] = V);
  return out[0];
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
  // sampler storage - average
  double m_w_avg;
  vec m_kappa_avg, m_inv_tau_sq_avg, m_beta_avg;
  mat m_inv_Sigma_avg;


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
    m_beta_draw.set_size(m_delta.n_elem);
    m_beta_draw.zeros();
    m_kappa_draw.set_size(m_delta.n_elem);
    m_kappa_draw.ones();
    m_inv_tau_sq_draw.set_size(m_delta.n_elem);
    m_inv_tau_sq_draw.fill(1 / pow(10, 6));
    m_inv_Sigma.set_size(m_y_mx.n_cols, m_y_mx.n_cols);
    m_inv_Sigma.eye();
  }
  void initialize_avg () {
    m_w_avg = 0;
    m_beta_avg.copy_size(m_beta_draw);
    m_beta_avg.zeros();
    m_kappa_avg.copy_size(m_kappa_draw);
    m_kappa_avg.zeros();
    m_inv_tau_sq_avg.copy_size(m_inv_tau_sq_draw);
    m_inv_tau_sq_avg.zeros();
    m_inv_Sigma_avg.copy_size(m_inv_Sigma);
    m_inv_Sigma_avg.zeros();
  }
  // update functions
  void beta_update(const mat &chol_U_Inv_Sigma, const vec &y_vec, const double &blocks) {
    mat UX = chol_U_Inv_Sigma * m_X_alpha;
    vec U_y_vec = chol_U_Inv_Sigma * y_vec;
    double q = ceil(m_delta.n_elem / blocks);
    for (uword j = 0; j < blocks; j++) {
      mat UX_j, R_mx_j, inv_B_N_j, B_N_j, UX_shed;
      vec b_N_j, beta_draw_shed, end_vec(2);
      // create selectors for block
      int start = j * q;
      end_vec(0) = (j + 1) * q;
      end_vec(1) = m_delta.n_elem;
      int end = end_vec.min() - 1;
      // block quantities and sample
      UX_j = UX.cols(start, end);
      R_mx_j = m_GRR_mx.submat(start, start, end, end);
      inv_B_N_j = R_mx_j + UX_j.t() * UX_j;
      B_N_j = inv_sympd(inv_B_N_j);
      UX_shed = UX;
      UX_shed.shed_cols(start, end);
      beta_draw_shed = m_beta_draw;
      beta_draw_shed.shed_rows(start, end);
      b_N_j = B_N_j * UX_j.t() * (U_y_vec - UX_shed * beta_draw_shed);
      // draw
      m_beta_draw.subvec(start, end) = mvnrnd(b_N_j, B_N_j);
    }
  }
  void kappa_update () {
    vec beta_sq, sd_term, w_1, w_2, v;
    beta_sq = pow(m_beta_draw, 2);
    sd_term = .5 * m_inv_tau_sq_draw % beta_sq;
    w_1 = exp(log(1 - m_w) - .5 * log(m_epsilon) - sd_term / m_epsilon);
    w_2 = exp(log(m_w) - sd_term);
    // draw
    v = r_binom(m_beta_draw.n_elem, w_1 / (w_1 + w_2));
    m_kappa_draw = m_epsilon * v + (1 - v);
  }
  void inv_tau_sq_update () {
  vec beta_sq, a_N_2;
  beta_sq = pow(m_beta_draw, 2);
  a_N_2 = m_a_2 + .5 * beta_sq / m_kappa_draw;
  // draw
  m_inv_tau_sq_draw = r_gamma(m_inv_tau_sq_draw.n_elem, m_a_1 + .5, a_N_2);
}
  void w_update() {
    double v, V;
    v = 1 + sum(m_kappa_draw == 1);
    V = 1 + sum(m_kappa_draw != 1);
    // draw
    m_w = r_beta(v, V);
  }
  void GRR_update() { m_GRR_mx = diagmat(m_inv_tau_sq_draw / m_kappa_draw); }
  void inv_Sigma_update (const vec &y_vec) {
    double v_N;
    vec r_vec;
    mat S_N, inv_S_N, r_mx;
    r_vec = y_vec - m_X_alpha * m_beta_draw;
    r_mx.set_size(m_X_block.n_rows, m_y_mx.n_cols);
    for (uword i = 0; i < m_y_mx.n_cols; i++) {
      r_mx.col(i) = r_vec.subvec(i * m_X_block.n_rows, (i + 1) * m_X_block.n_rows - 1);
    }
    v_N = m_v_0 + m_X_block.n_rows + m_y_mx.n_cols;
    S_N = m_S_0 + r_mx.t() * r_mx;
    inv_S_N = inv_sympd(S_N);
    // draw
    m_inv_Sigma = wishrnd(inv_S_N, v_N);
  }
  // STORE SAMPLES
  void update_avg() {
    m_beta_avg += m_beta_draw;
    m_kappa_avg += m_kappa_draw;
    m_inv_tau_sq_avg += m_inv_tau_sq_draw;
    m_w_avg += m_w;
    m_inv_Sigma_avg += m_inv_Sigma;
  }
  void normalize_avg(const int &iterations, const int &burn) {
    double n = iterations - burn;
    m_beta_avg /= n;
    m_kappa_avg /= n;
    m_inv_tau_sq_avg /= n;
    m_w_avg /= n;
    m_inv_Sigma_avg /= n;
  }


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
  void sample(const bool &RATS, const double &blocks, const int &iterations, const int &burn, const bool &point) {
    List out(5);
    // helpers
    mat Kron_hlpr(m_X_block.n_rows, m_X_block.n_rows, fill::eye);
    // hyperparameters Wishart
    set_Wishart_hp(RATS);
    // initial draws
    set_initial_draws();
    initialize_avg();
    GRR_update();
    mat chol_U_inv_Sigma = chol(m_inv_Sigma);
    mat chol_U_Inv_Sigma = kron(chol_U_inv_Sigma, Kron_hlpr);
    vec y_vec = vectorise(m_y_mx);
    for (size_t iter = 0; iter < iterations; iter++) {
      beta_update(chol_U_Inv_Sigma, y_vec, blocks);
      kappa_update();
      inv_tau_sq_update();
      GRR_update();
      w_update();
      inv_Sigma_update(y_vec);
      chol_U_inv_Sigma = chol(m_inv_Sigma);
      chol_U_Inv_Sigma = kron(chol_U_inv_Sigma, Kron_hlpr);
      if (point) {
        if (iter >= burn) {
          update_avg();
        }
      }
    }
    if (point) {
      normalize_avg(iterations, burn);
    }
  }
  List prep_output (const bool point) {
    List out;
    if (point) {
      out = List::create(
        _["beta"] = m_beta_avg,
        _["kappa"] = m_kappa_avg,
        _["tau.sq"] = 1 / m_inv_tau_sq_avg,
        _["omega"] = m_w_avg,
        _["Sigma"] = inv_sympd(m_inv_Sigma_avg)
      );
    } else {

    }
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
List JBU_sample(const mat &X_block, const mat &y_mx, const int &K,
  const bool RATS = true, const double blocks = 1,
  const int iterations = 1000, const int burn = 200,
  const bool point = true
) {
  List out;
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  mod.fill_X_alpha();
  mod.sample(RATS, blocks, iterations, burn, point);
  out = mod.prep_output(point);
  return out;
}
