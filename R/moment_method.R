moment_vanila <- function(X, Y){
  n = length(Y)
  XtrY = crossprod(X, Y)
  V_hat = t(X) %*% X
  eigen.V = eigen(V_hat)
  p = nrow(V_hat)
  rank_r = sum(eigen.V$values>1e-8)
  ii = which(eigen.V$values>1e-8)
  W =   eigen.V$vectors %*% diag( c( (eigen.V$values[ii])^(-1), rep(0,p - rank_r)   )  ) %*% t(eigen.V$vectors)
  h2 =  t(XtrY) %*% W %*% XtrY/ sum(Y^2)  - rank_r/n
  return(c(h2, (n/(n - rank_r))^2 * ((2* rank_r ) * ( 1- (h2)^2 )/ n  + 4*(h2)^2 ) * ( 1- h2 )/n ))
  
}
