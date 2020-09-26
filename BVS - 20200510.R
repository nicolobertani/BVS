# Bayesian Variable Selection

# Setup -------------------------------------------------------------------

rm(list = ls())
if (!requireNamespace("Matrix")) install.packages('Matrix')
library(Matrix)
if (!requireNamespace("MASS")) install.packages('MASS')
library(MASS)
# if (!requireNamespace("lars")) install.packages('lars')
# library(lars)



# Simulate data ---------------------------------------------------------

n.regressors <- 12
train.size <- 1e3
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
  
# plot
par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = c(2,2,0,0))
for(i in seq(n.regressors)) {
  plot(y ~ X.std[, i], xaxt = 'n', yaxt = 'n')
  mtext(paste0(expression(beta), " = ", round(beta.true[i], 2)), 1, .5)
}



# Discrete Bayesian model selection ---------------------------------------

# define the function for model evidence
model.evidence <- function(y.model, X.model, a_0, A_0, b_0, c_0) {
  # create parameters
  n <- length(y.model)
  inv.A_0 <- chol2inv(chol((A_0)))
  inv.A <- t(X.model) %*% X.model + chol2inv(chol(A_0))
  A <- chol2inv(chol(inv.A))
  a <- A %*% (t(X.model) %*% y.model + inv.A_0 %*% a_0)
  b <- b_0 + n / 2
  c <- c_0 + (t(y.model) %*% y.model - t(a) %*% inv.A %*% a + t(a_0) %*% inv.A_0 %*% a_0) / 2
  
  log.evidence <- - n / 2 * log(2 * pi) + 
    log(det(A)) / 2 - log(det(A_0)) / 2 +
    b_0 * log(c_0) - b * log(c) +
    lgamma(b) - lgamma(b_0)
  return(as.numeric(log.evidence))
}

# create list of all possible models
all.possible.models <- as.matrix(expand.grid(replicate(n.regressors, list(0,1), simplify = F)))[-1, ]
all.possible.models.list <- lapply(seq(nrow(all.possible.models)), function(i) {unlist(all.possible.models[i, ])})

# compute model evidence for all possible models
evidence.results <- unlist(lapply(all.possible.models.list, function(ar) {
  model.evidence(y.train, X.train.std[, as.logical(ar)], rep(0, sum(ar)), diag(1e6, sum(ar)), .001, .001)
}))

# model with highest evidence
max(evidence.results)
which.max(evidence.results)
selected.m <- all.possible.models.list[[which.max(evidence.results)]]
rbind("truth" = beta.true != 0,
      "evidence" = all.possible.models.list[[which.max(evidence.results)]])


# estimate the model, available in closed-form
model.estimates <- function(y.model, X.model, a_0, A_0, b_0, c_0) {
  n <- length(y.model)
  inv.A_0 <<- chol2inv(chol((A_0)))
  inv.A <<- t(X.model) %*% X.model + inv.A_0
  A <- chol2inv(chol(inv.A))
  a <<- A %*% (t(X.model) %*% y.model + inv.A_0 %*% a_0)
  b <- b_0 + n / 2
  c <- c_0 + (t(y.model) %*% y.model - t(a) %*% inv.A %*% a + t(a_0) %*% inv.A_0 %*% a_0) / 2
  return(list(
    a = a,
    A = A,
    b = b, 
    c = c
  ))
}
selected.m <- append(list(active.beta = selected.m),
  model.estimates(y.train, X.train.std[, as.logical(selected.m)], rep(0, sum(selected.m)), diag(1e6, sum(selected.m)), .001, .001)
)
BMS.beta <- rep(0, n.regressors)
BMS.beta[as.logical(selected.m$active.beta)] <- selected.m$a
rbind(BMS.beta, beta.true)

par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = rep(.2, 4))
for(i in seq(n.regressors)) {
  plot(y ~ X.std[, i], xaxt = 'n', yaxt = 'n')
  abline(a = 0, b = beta.true[i], col = 2, lwd = 2)
  abline(a = 0, b = BMS.beta[i], col = 3, lwd = 2)
  if(i == 1) {
    legend("topleft", legend = c("Truth", "BMS"), col = c(2,3), lwd = 2)
  }
}



# Bayesian Ridge ----------------------------------------------------------

# estimation require Gibbs sampler
# define update functions
update.beta <- function(y.model, X.model, sigma.sq, psi) {
  # posterior parameters for regressors
  inv.A_0 <- diag(ncol(X.model)) / psi
  inv.A <- t(X.model) %*% X.model + inv.A_0
  A <- chol2inv(chol(inv.A))
  a <- A %*% (t(X.model) %*% y.model)
  # draw
  output <- mvrnorm(1, a, sigma.sq * A)
  return(output)
}
update.sigma.sq <- function(y.model, X.model, beta, b_0, c_0) {
  n <- length(y.model)
  inv.sigma.sq <- rgamma(1,
                         b <- b_0 + n / 2,
                         c <- c_0 + t(y.model - X.model %*% beta) %*% (y.model - X.model %*% beta) / 2
  )
  return(1 / inv.sigma.sq)
}
update.psi <- function(beta, sigma.sq, s_0, t_0) {
  inv.psi <- rgamma(1,
                    s <- s_0 + n.regressors / 2,
                    t <- t_0 + t(beta) %*% beta / (2 * sigma.sq)
  )
  return(1 / inv.psi)
}

# prepare sampler
n.iter <- 700
burn.in.ratio <- 2 / 7
beta.draws <- matrix(NA, ncol = n.regressors, nrow = n.iter * (1 - burn.in.ratio))


# Run Gibbs sampler
iteration <- 0
saving.iteration <- 0
sigma.sq.draw <- 1
psi.draw <- 1

while (iteration < n.iter) {
  iteration <- iteration + 1
  if(iteration == 1) pb <- txtProgressBar(min = 0, max = n.iter, initial = 0, style = 3)
  
  beta.draw <- update.beta(y.train, X.train.std, sigma.sq.draw, psi.draw)
  sigma.sq.draw <- update.sigma.sq(y.train, X.train.std, beta.draw, .001, .001)
  psi.draw <- update.psi(beta.draw, sigma.sq.draw, 1e1, 1e1)
  
  if (iteration > (n.iter * burn.in.ratio)) {
    saving.iteration <- saving.iteration + 1 
    beta.draws[saving.iteration, ] <- beta.draw
    }
  
  setTxtProgressBar(pb, iteration)
}

# plot
par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = c(2,2,0,0))
for(i in seq(n.regressors)) {
  plot(beta.draws[, i], xaxt = 'n', yaxt = 'n', type = 'l')
}

BR.beta <- colMeans(beta.draws)
rbind( 
  "Bayesian Ridge" =  BR.beta,
  "Model selection" = BMS.beta,
  "True model" = beta.true
)

# plot
par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = rep(.2, 4))
for(i in seq(n.regressors)) {
  plot(y ~ X.std[, i], xaxt = 'n', yaxt = 'n')
  abline(a = 0, b = beta.true[i], col = 2, lwd = 2)
  abline(a = 0, b = BMS.beta[i], col = 3, lwd = 2)
  abline(a = 0, b = BR.beta[i], col = 4, lwd = 2)
  if(i == 1) {
    legend("topleft", legend = c("Truth", "BMS", "Bayesian Ridge"), col = c(2,3,4), lwd = 2)
  }
}




# Dirac Spike and Slab ----------------------------------------------------

# estimation requires Gibbs samples
# define update functions
ratio <- function(delta.index, delta, sigma.sq, psi) {
  d.with <- delta
  d.with[delta.index] <- 1
  d.wout <- delta
  d.wout[delta.index] <- 0
  X.with <- X.train.std[, as.logical(d.with)]
  X.wout <- X.train.std[, as.logical(d.wout)]
  n <- length(y.train)
  
  # parameters with
  k.with <- sum(d.with)
  inv.A.N.with <- t(X.with) %*% X.with + 1 / psi * diag(k.with)
  A.N.with <- solve(as.matrix(inv.A.N.with))
  a.with <- A.N.with %*% t(X.with) %*% y.train
  
  # manage case of all deltas being zero
  if (sum(d.wout) == 0) {
    
    output <- - (t(a.with) %*% inv.A.N.with %*% a.with) / sigma.sq + log(psi / det(A.N.with))
    
  } else {
    
    # parameters without
    k.wout <- sum(d.wout)
    inv.A.N.wout <- t(X.wout) %*% X.wout + 1 / psi * diag(sum(d.wout))
    A.N.wout <- solve(as.matrix(inv.A.N.wout))
    a.wout <- A.N.wout %*% t(X.wout) %*% y.train
    output <- - (t(a.with) %*% inv.A.N.with %*% a.with - t(a.wout) %*% inv.A.N.wout %*% a.wout) / sigma.sq + log(psi * det(A.N.wout) / det(A.N.with))
  }
  
  ratio.out <- exp(output / 2)
  return(ratio.out)
}


n.iter <- 2500
burn.in.ratio <- 1/5
draws <- matrix(NA, ncol = 2 * n.regressors + 2, nrow = n.iter * (1 - burn.in.ratio))

n <- length(y.train)
# hyperparameters
b_0 <- .001
c_0 <- .001
s_0 <- 1
t_0 <- 1
# omega.draw <- rbeta(1, s_0, t_0)
omega.draw <- .5
psi.draw <- 1e4 # inverse shrinkage parameter
sigma.sq.draw <- 1

# draw non-zero regressors
alpha <- rbinom(n.regressors, 1, omega.draw)
{
  iteration <- 0
  saving.iteration <- 0
  while (iteration < n.iter) {
    if(iteration == 0) pb <- txtProgressBar(min = 0, n.iter, label = 'Gibbs Sampler', style = 3)
    iteration <- iteration + 1
    
    alpha.positions <- sample(seq(n.regressors))
    for(i in alpha.positions) {
      alpha[i] <- rbinom(1, 1, 1 / (1 + ratio(i, alpha, sigma.sq.draw, psi.draw) * (1 - omega.draw) / omega.draw))
    }
    # all.equal(
    #   sapply(alpha.positions, function(i) {ratio(i, alpha, sigma.sq.draw, psi.draw)}),
    #   sapply(alpha.positions, function(i) {rcpp_ratio(alpha, i, X.train.std, y.train, sigma.sq.draw, psi.draw)})
    # )
    # microbenchmark::microbenchmark(
    #   sapply(alpha.positions, function(i) {ratio(i, alpha, sigma.sq.draw, psi.draw)}),
    #   sapply(alpha.positions, function(i) {rcpp_ratio(alpha, i, X.train.std, y.train, sigma.sq.draw, psi.draw)})
    # )
    
    if (sum(alpha) != 0) {
      # normal distribution of regressors
      X.alpha <- X.train.std[, as.logical(alpha)]
      inv.A.alpha <- t(X.alpha) %*% X.alpha + 1 / psi.draw * diag(sum(alpha))
      A.alpha <- solve(as.matrix(inv.A.alpha))
      a.alpha <- A.alpha %*% (t(X.alpha) %*% y.train)
      # sample sigma sq
      b <- b_0 + (n - 1) / 2
      c <- c_0 + as.numeric(t(y.train) %*% y.train - t(a.alpha) %*% inv.A.alpha %*% a.alpha) / 2
      sigma.sq.draw <- 1 / rgamma(1, shape = b, rate = c)
      # draw active regressors
      beta.alpha.draw <- mvrnorm(1, a.alpha, sigma.sq.draw * A.alpha)
      # draw regressors
      beta.draw <- rep(0, n.regressors)
      beta.draw[as.logical(alpha)] <- beta.alpha.draw
      
    } else {
      
      # sample sigma sq
      b <- b_0 + (length(y.train) - 1) / 2
      c <- c_0 + as.numeric(t(y.train) %*% y.train) / 2
      sigma.sq.draw <- 1 / rgamma(1, shape = b, rate = c)
      # draw regressors
      beta.draw <- rep(0, n.lags)
    }
    
    omega.draw <- rbeta(1, s_0 + sum(alpha), t_0 - sum(alpha) + length(alpha))
    
    # save draws
    if (iteration > (n.iter * burn.in.ratio)) {
      saving.iteration <- saving.iteration + 1
      draws[saving.iteration, ] <- c(alpha, beta.draw, sigma.sq.draw, omega.draw)
    }
    setTxtProgressBar(pb, iteration)
  }
}
Rcpp::sourceCpp("~/OneDrive - INSEAD/Github/BVS code/Dirac_SS_IP.cpp")
output <- rcpp_Dirac_SS(X.train.std, y.train, 2e3, 5e2, update_psi = 0, fixed_psi = length(y.train))
round(rbind(apply(draws[, 1:12], 2, mean),
            apply(output[[1]], 2, mean)), 2)
round(rbind(apply(draws[, 1:12 + 12], 2, mean),
            apply(output[[2]], 2, mean)), 2)
c(mean(draws[, 25]), mean(output[[3]]))
mean(output[[4]])
c(mean(draws[, 26]), mean(output[[5]]))
Rcpp::sourceCpp("~/OneDrive - INSEAD/Github/BVS code/Dirac_SS_gP.cpp")
output_gP <- rcpp_Dirac_SS_g(X.train.std, y.train, 2e3, 5e2, fixed_g = length(y.train))
round(rbind(apply(output[[1]], 2, mean),
            apply(output_gP[[1]], 2, mean)), 2)
round(rbind(apply(output[[2]], 2, mean),
            apply(output_gP[[2]], 2, mean)), 2)
round(rbind(apply(output[[3]], 2, mean),
            apply(output_gP[[3]], 2, mean)), 2)
c(mean(output[[4]]), mean(output_gP[[4]]))
c(mean(output[[5]]), mean(output_gP[[5]]))


# active barplot
par(mfrow = c(1,1), oma = rep(0,4) + .1, mar = rep(.2, 4))
barplot(colMeans(draws[, seq(n.regressors)]), col = (beta.true != 0) + 2, border = NA)

# plot
par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = rep(.2, 4))
for(i in seq(n.regressors)) {
  plot(y ~ X.std[, i], xaxt = 'n', yaxt = 'n')
  abline(a = 0, b = beta.true[i], col = 2, lwd = 2)
  abline(a = 0, b = BMS.beta[i], col = 3, lwd = 2)
  abline(a = 0, b = BR.beta[i], col = 4, lwd = 2)
  abline(a = 0, b = colMeans(draws[, seq(n.regressors) + n.regressors])[i], col = 5, lwd = 2)
  if(i == 1) {
    legend("topleft", legend = c("Truth", "BMS", "Bayesian Ridge", "Dirac"), col = c(2,3,4, 5), lwd = 2)
  }
}



# ACSS - NMIG  ------------------------------------------------------------

update.gamma <- function() {
  # posterior parameters for mixture
  beta.sq <- as.vector(beta.draw ^ 2)
  w_1 <- exp(log(1 - w) - log(epsilon) / 2 - beta.sq / (2 * epsilon * tau.sq))
  w_2 <- exp(log(w) - beta.sq / (2 * tau.sq))
  # draw
  v <- rbinom(n.regressors, 1, as.vector(w_1 / (w_1 + w_2)))
  output <- v * epsilon + (1 - v)
  return(output)
}

update.inv.tau.sq <- function() {
  beta.sq <- as.vector(beta.draw ^ 2)
  inv.tau.sq <- rgamma(n.regressors, 
                       a_1 + .5,
                       a_2 + beta.sq / (2 * kappa))
  return(inv.tau.sq)
}

update.omega <- function() {
  w <- rbeta(1,
             1 + sum(kappa == 1),
             1 + sum(kappa == epsilon))
  return(w)
}

update.inv.A.matrix <- function() {
  diag(inv.tau.sq / kappa)
}

update.beta <- function() {
  # posterior parameters for regressors
  inv.B.N <- as.matrix(inv.A_0 + t(X.train.std) %*% X.train.std * inv.sigma.sq)
  # B.N <- solve(inv.B.N) # solve
  t.chol.B.N <- solve(chol(inv.B.N))
  B.N <- t.chol.B.N %*% t(t.chol.B.N)
  b.N <- B.N %*% t(X.train.std) %*% y.train * inv.sigma.sq
  # draw
  output <- mvrnorm(1, b.N, B.N)
  return(output)
}

update.inv.sigma.sq <- function() {
  inv.sigma.sq <- rgamma(1,
                         b_1 + n / 2,
                         b_2 + norm(y.train - X.train.std %*% beta.draw, '2') ^ 2 / 2
  )
  return(inv.sigma.sq)
}

# storage
beta.draws.ACSS <- matrix(NA, ncol = n.regressors, nrow = n.iter * (1 - burn.in.ratio))
gamma.draws.ACSS <- matrix(NA, ncol = n.regressors, nrow = n.iter * (1 - burn.in.ratio))

# Hyperparameters
epsilon <- .001 # small value 
a_1 <- 5 # shape of tau prior
a_2 <- 50 # rate of tau prior
w <- 1 / 2 # indifference prior for probability of shrinkage
b_1 <- 1e-4 # shape of prior for inverse sigma sq
b_2 <- 1e-4 # rate of prior for inverse sigma sq

# Initial values for sampler
kappa <- rep(1, n.regressors)
tau.sq <- rep(1e6, n.regressors)
inv.sigma.sq <- 1
inv.A_0 <- diag(1 / (kappa * tau.sq))
beta.draw.active <- rep(0, n.regressors)

# Run Gibbs sampler
iteration <- 0
saving.iteration <- 0
while (iteration < n.iter) {
  iteration <- iteration + 1
  if(iteration == 1) pb <- txtProgressBar(min = 0, max = n.iter, initial = 0, style = 3)
  
  # draw beta
  beta.draw <- update.beta()
  # draw ridge matrix
  kappa <- update.gamma()
  inv.tau.sq <- update.inv.tau.sq()
  tau.sq <- 1 / inv.tau.sq
  inv.A_0 <- update.inv.A.matrix()
  # draw w
  w <- update.omega()
  # draw sigma.sq
  inv.sigma.sq <- update.inv.sigma.sq()
  
  if (iteration > (n.iter * burn.in.ratio)) {
    saving.iteration <- saving.iteration + 1    
    beta.draws.ACSS[saving.iteration, ] <- as.vector(beta.draw)
    gamma.draws.ACSS[saving.iteration, ] <- as.vector(kappa)
  }
  
  setTxtProgressBar(pb, iteration)
}


# active barplot
par(mfrow = c(1,1), oma = rep(0,4) + .1, mar = rep(.2, 4))
barplot(colMeans(gamma.draws.ACSS == 1), col = (beta.true != 0) + 2, border = NA)

# plot
par(mfrow = c(3,4), oma = rep(0,4) + .1, mar = rep(.2, 4))
for(i in seq(n.regressors)) {
  plot(y ~ X.std[, i], xaxt = 'n', yaxt = 'n')
  abline(a = 0, b = beta.true[i], col = 2, lwd = 2)
  abline(a = 0, b = colMeans(draws[, seq(n.regressors) + n.regressors])[i], col = 5, lwd = 2)
  abline(a = 0, b = colMeans(beta.draws.ACSS)[i], col = 'pink', lwd = 2)
  if(i == 1) {
    legend("topleft", legend = c("Truth", "Dirac", "ACSS"), col = c(2, 5, "pink"), lwd = 2)
  }
}
































