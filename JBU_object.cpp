#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


class JBU_model {

  mat m_X_block;
  mat m_y_mx;
  double m_psi = pow(10, 6);

public:
  // construct and print X and y
  void set_X_block (mat value) { m_X_block = value ;}
  mat get_X_block () { return m_X_block; }
  void set_y_mx (mat value) { m_y_mx = value ;}
  mat get_y_mx () { return m_y_mx; }

  // approximate Dirac filter
  void filter () {
    for (size_t i = 0; i < m_y_mx.n_cols; i++) {
      for (size_t j = 0; j < m_X_block.n_cols; j++) {
        double inv_A_N = dot(m_X_block.col(j), m_X_block.col(j)) + 1 / m_psi;
        double A_N = 1 / inv_A_N;
        double a_N = A_N * dot(m_X_block.col(j), m_y_mx.col(i));
        double B = dot(m_y_mx.col(i), m_y_mx.col(i)) - pow(a_N, 2) * inv_A_N;
        cout << log(A_N) - (m_X_block.n_rows - 1) * log(B) << "\n";
      }
      /* code */
    }
  }

};

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
void test_assignments (const mat &X_block, const mat &y_mx) {
  JBU_model mod;
  mod.set_X_block(X_block);
  mod.set_y_mx(y_mx);
  mod.filter();
}
