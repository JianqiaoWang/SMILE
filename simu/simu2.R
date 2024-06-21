library(dplyr)
n = 800
p = 1000
h2 = 0.4

devtools::load_all("/n/holystore01/LABS/xlin/Lab/wangjq/SMILE")
res = list()
sim = function(i){
  set.seed(i)
  print(i)
  X = matrix(rnorm(n*p), n, p) %*% Sigma_sqrt
  S <-  X %*% t(X)
  S = (n/sum(diag(S) )) * S
  S_eigen = eigen(S, only.values = T)
  return(S_eigen$values)
}

for(rho in ((1:9)/10)){
Sigma = rho^(abs(outer(1:p, 1:p, "-")))
Sigma_sqrt = with(eigen(Sigma), vectors %*% (values^(1/2) * t(vectors)))
res[[as.character(rho) ]] = lapply(1:500, sim )
}

saveRDS(res,  file = "/n/holystore01/LABS/xlin/Lab/wangjq/SMILE/simu/eigenvalue.rds")
