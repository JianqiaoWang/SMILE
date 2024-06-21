
library(dplyr)
library(parallel)
p <- 1000
h2 = 0.4

n <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
if(is.na(n)){n = 800}

devtools::load_all("/n/holystore01/LABS/xlin/Lab/wangjq/SMILE")
ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD <- readRDS("/n/holystore01/LABS/xlin/Lab/wangjq/svd_heri/prepare/ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD.rds")
Sigma = ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD[1:1000, 1:1000]
eig1 = eigen(Sigma)
Sigma_sqrt = eig1$vectors %*% (eig1$values^(1/2) * t(eig1$vectors))

beta_coef = c( eig1$vectors[,1])
beta_coef[-sample(1:1000, 20)] = 0
beta_coef = beta_coef/sqrt( as.numeric( t( beta_coef) %*% Sigma %*% beta_coef ) ) * sqrt(h2)
blist = list( c(1:(p/2) ) , c((p/2) + 1:(p/2) ) ) # define a blist

W1 = NULL # diagonal
W2= solve(Sigma) # true inverse
#W5 = diag( as.vector(scale(1/c(eig1$vectors[,1]))))

n_E = n
X_E = matrix(rnorm(n_E *p), n_E, p) %*% Sigma_sqrt
# specify the inverse
out_of_sample_inverse =  as.matrix(Matrix::bdiag( mat_inv( cov(X_E[, 1:500]) ),
                                                  mat_inv( cov(X_E[, 501:1000]) )  ))
W4 = out_of_sample_inverse


 mat_inv = function(mat){

   with(eigen( mat), {

        significant_indices <- which(values > 1e-8)

        vectors[, significant_indices] %*% diag(1 / values[significant_indices]) %*% t(vectors[, significant_indices])
   }
        )

 }

# W6 =  mat_inv(cov(X_E))




sim = function(i){

  set.seed(i)
  print(i)

  X = matrix(rnorm(n*p), n, p) %*% Sigma_sqrt

  # generate the outcome
  eps = rnorm(n) * sqrt(1 - h2)
  y = X %*% beta_coef + eps

  # specify the inverse
  in_sample_inverse =  as.matrix(Matrix::bdiag( mat_inv( cov(X[, 1:500]) ),
                                                mat_inv( cov(X[, 501:1000]) )  ))


  W3= in_sample_inverse

  # calculate the heritability
  res1 = scoring_vanila(X, y, W = W1, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)
  res2 = scoring_vanila(X, y, W = W2, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)
  res3 = scoring_vanila(X, y, W = W3, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)
  res4 = scoring_vanila(X, y, W = W4, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)

  #res5 = scoring_vanila(X, y, W = W5, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)

  return( list(res1, res2, res3, res4) )
}


cmpr <- mclapply( 1:500, sim , mc.cores = 5)

filename = paste0("/n/holystore01/LABS/xlin/Lab/wangjq/SMILE/simu/W_efficiency_", n, ".rds")

saveRDS(cmpr, file = filename)

# module load gcc/12.2.0-fasrc01
# source /n/home09/wangjq/intel/oneapi/mkl/2022.2.1/env/vars.sh intel64
# module load R
# Rfile='/n/holystore01/LABS/xlin/Lab/wangjq/SMILE/simu/simu3.R'
# sbatch --mem=50G -c 6 -t 0-12:10 -p xlin --wrap="/n/home09/wangjq/software/R-4.2.0/bin/Rscript $Rfile 800"
