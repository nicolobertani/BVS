#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

vec r_stdnorm(const int &n) {
  Function r_norm("rnorm");
  NumericVector tmp = r_norm(_["n"] = n);
  vec out(n);
  for (uword i = 0; i < n; i++) {
    out(i) = tmp[i];
  }
  return out;
}

// [[Rcpp::export]]
vec rmvn_eig(vec &mu, mat &Sigma) {
  mat U, D_sqrt;
  vec d, draw, out;
  eig_sym(d, U, Sigma);
  D_sqrt = diagmat(sqrt(d));
  draw = r_stdnorm(mu.n_elem);
  out = mu + U * D_sqrt * draw;
  return out;
}

// [[Rcpp::export]]
vec rmvn_chol(vec &mu, mat &Sigma) {
  return mvnrnd(mu, Sigma);
}

// [[Rcpp::export]]
vec rcpp_mvrnorm(const vec &mu, const mat &Sigma) {
    vec out(mu.n_elem);
    Environment pkg = Environment::namespace_env("MASS");
    Function f = pkg["mvrnorm"];
    NumericVector tmp = f(1, mu, Sigma);
    for (uword i = 0; i < mu.n_elem; i++) {
      out(i) = tmp[i];
    }
    return out;
}

// [[Rcpp::export]]
mat rwish (int v, mat S) {
  return wishrnd(S, v);
}
