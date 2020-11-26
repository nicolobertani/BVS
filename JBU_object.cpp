#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


class JBU_model {

  mat m_X_block;
  mat m_y_mx;
  double m_psi = pow(10, 4);
  umat m_delta;

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

public:
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
      m_delta.col(i) = sort(order_v.head(K));
    }
  }
  // approximate Dirac filter: return m_delta
  umat get_delta_mx () { return m_delta; }
};

/*
END OF JBU MODEL CLASS
START R FUNCTIONS
*/

// [[Rcpp::export]]
mat test_function_X (const mat &X_block, const mat &y_mx) {
  JBU_model mod;
  mod.set_X_block(X_block);
  return mod.get_X_block();
}

// [[Rcpp::export]]
mat test_function_y (const mat &X_block, const mat &y_mx) {
  JBU_model mod;
  mod.set_y_mx(y_mx);
  return mod.get_y_mx();
}

// [[Rcpp::export]]
umat JBU_filter (const mat &X_block, const mat &y_mx, const int &K) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter(K);
  return mod.get_delta_mx() + 1;
}
