library(RcppArmadillo)
setwd("~/OneDrive - INSEAD/Github/BVS code/")
Rcpp::sourceCpp("Dirac_active_coefficients_sampler.cpp")



# simulate data -----------------------------------------------------------

# Simulate data ---------------------------------------------------------

n.regressors <- 12
train.size <- 1e2
test.size <- 1e2
data.size <- train.size + test.size

# simulate covariates
X <- matrix(rnorm(n.regressors * data.size), ncol = n.regressors)
dim(X)
X.train <- head(X, train.size)
X.test <- tail(X, test.size)

# standardize
X.m <- apply(X.train, 2, mean)
X.sd <- apply(X.train, 2, sd)
# # checks
# apply((t(X.train) - X.m) / X.sd, 1, sum)
# apply((t(X.train) - X.m) / X.sd, 1, sd)

X.std <- t((t(X) - X.m) / X.sd)
X.train.std <- head(X.std, train.size)
X.test.std <- tail(X.std, test.size)
# # checks
# apply(X.train.std, 2, sum)
# apply(X.train.std, 2, sd)

# simulate beta
# we will assume that 6 out of the 12 possible regressors are active
beta.true <- runif(n.regressors, -1, 1)
beta.true <- beta.true * sample(c(rep(0, n.regressors / 2), rep(1, n.regressors / 2)))

# simulate the observations
# assuming standard normal noise
y <- X.std %*% beta.true + rnorm(data.size)
y.train <- head(y, train.size)
y.test <- tail(y, test.size)


# run ---------------------------------------------------------------------

delta <- rbinom(n.regressors, 1, .5)
# rcpp_ratio(delta, 11, mx, 1:10, 2, .5)
output <- rcpp_Dirac_SS(X.train.std, y.train, 10000, 2000)
apply(output[[1]], 2, mean)

