#ifndef _SAMPLING_FNS
#define _SAMPLING_FNS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double r_beta(const double &v, const double &V) {
  Rcpp::Function r_beta("rbeta");
  Rcpp::NumericVector out = r_beta(1, Rcpp::_["shape1"] = v, Rcpp::_["shape2"] = V);
  return out[0];
}

#endif
