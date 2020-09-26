## High performance BVS

The software in this repository efficiently implements a subset of the Bayesian Variable Selection (BVS) techniques generally referred to as Spike-and-Slab regression.

More specifically:

* Dirac_SS_IP implements Dirac Spike-and-Slab regression with the independent prior.

* Dirac_SS_gP implements Dirac Spike-and-Slab regression with the g-prior.

Estimation relies on Gibbs sampling.

For a review of the topic, refer to:
Malsiner-Walli, G. and Wagner, H. (2018). Comparing spike and slab priors for bayesian variable selection. arXiv preprint [arXiv:1812.07259](https://arxiv.org/abs/1812.07259).

The software is written in C++ for performance, and is meant for usage in R via [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html).

This project is ongoing. Contributions are welcome.

### License

High performance Bayesian Variable Selection for R using C++ via Rcpp and RcppArmadillo

Copyright (C) 2020  Nicolò Bertani

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>. You should also find a copy of this license in this repository.

Please use Github Issues to get in touch.
