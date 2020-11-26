#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


class JBU_model {

  mat m_X_block, m_y_mx;
  umat m_delta;
  sp_mat m_X_alpha;
  // hyperparameters
  double m_psi = pow(10, 4); // beta
  double a_1 = 5; // tau
  double a_2 = 50; // tau
  double epsilon = .005; // gamma
  double m_w = .5; // omega
  double m_v_0;
  mat m_S_0;

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
  void sample(const bool &RATS, const bool &large) {
    // hyperparameters Wishart
    set_Wishart_hp(RATS);
    cout << m_v_0 << "\n";
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
void test(const mat &X_block, const mat &y_mx, const int &K, const bool &RATS, const bool &large) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  mod.fill_X_alpha();
  mod.sample(RATS, large);
}
