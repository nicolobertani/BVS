#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec approx_pip(const mat &X_block, const mat &y_mx, const double psi, const int n_lags, const int n_bottom_series) {
  int P = n_lags * pow(n_bottom_series, 2);
  vec res(P);
  for (size_t i = 0; i < P; i++) {
    /* code */
  }
  return 0;
}
