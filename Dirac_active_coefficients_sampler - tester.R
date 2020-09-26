library(RcppArmadillo)
setwd("~/OneDrive - INSEAD/Github/BVS code/")
Rcpp::sourceCpp("Dirac_active_coefficients_sampler.cpp")


# run ---------------------------------------------------------------------

delta <- rbinom(12, 1, .5)
mx <- matrix(rep(1:12, each = 10), ncol = 12)
rcpp_ratio(delta, 12, mx, 1:10, 2, .5)

