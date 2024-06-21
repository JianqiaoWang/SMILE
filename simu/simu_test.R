library(dplyr)
n = 800
p = 1000
h2 = 0.4

beta_coef = rnorm(p)

devtools::load_all()

res = list()

Sigma = 0.8^(abs(outer(1:p, 1:p, "-")))

X = mvtnorm::rmvnorm(n, sigma = Sigma)
beta_coef = beta_coef/sqrt(t( beta_coef) %*% Sigma %*% beta_coef ) * sqrt(h2)

# generate the outcome
eps = rnorm(n) * sqrt(1 - h2)
y = X %*% beta_coef + eps

# calculate the heritability
res = scoring_vanila(X, y, W = NULL, theta0 =c(0.5, 0.5), max_iter = 100, eps = 1e-4, verbose = F)


W = diag(p)
Trace_XWX = if(is.null(W)){
  sum(X^2)
}else{
  sum( W * (t(X) %*% X)  )
}
W= n* W/Trace_XWX
res2 = scoring_conjugate(X, as.vector(y), kernel_fun = kernal_base, W = W, theta0 =c(0.5, 0.5), max_iter = 100, eps = 1e-4, verbose = F)

res3 = scoring_conjugate_AI(X, as.vector(y), kernel_fun = kernal_base, W = W, theta0 =c(0.5, 0.5), max_iter = 100, eps = 1e-4, verbose = T)


library(bigstatsr)
library(bigsnpr)
X = as_FBM(X)
res3 = scoring_conjugate_geno(X, as.vector(y), kernel_fun = kernal_FBM, W = W, theta0 =c(0.5, 0.5), max_iter = 100, eps = 1e-4, verbose = F)
