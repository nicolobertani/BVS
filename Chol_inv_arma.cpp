#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
void test () {
  mat A(5, 5, fill::randu);
  mat B = symmatu(A);
  mat C = symmatl(A);
  B.diag().fill(5);
  C.diag().fill(5);
  mat cholB = chol(B);
  mat cholC = chol(C);
  mat invcholB = inv(trimatu(cholB));
  mat invcholC = inv(trimatu(cholC));
  std::cout << approx_equal(inv_sympd(B), invcholB * invcholB.t(), "absdiff", .001) << "\n";
  std::cout << approx_equal(inv_sympd(C), invcholC * invcholC.t(), "absdiff", .001) << "\n";
}
