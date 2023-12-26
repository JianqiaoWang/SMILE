library(dplyr)
n = 800
p = 1000
h2 = 0.4

beta_coef = rnorm(p)

devtools::load_all()

res = list()

sim = function(i){

  X = mvtnorm::rmvnorm(n, sigma = Sigma)
  beta_coef = beta_coef/sqrt(t( beta_coef) %*% Sigma %*% beta_coef ) * sqrt(h2)

  # generate the outcome
  eps = rnorm(n) * sqrt(1 - h2)
  y = X %*% beta_coef + eps

  # calculate the heritability
  res = scoring_vanila(X, y, W = NULL, theta0 =c(0.5, 0.5), max_iter = 100, eps = 1e-4, verbose = F)

  return(res)
}

for(rho in ((1:9)/10)){
Sigma = rho^(abs(outer(1:p, 1:p, "-")))
res[[as.character(rho) ]] = lapply(1:100, sim )
}
