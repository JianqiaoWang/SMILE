
# SMILE

<!-- badges: start -->
<!-- badges: end -->

The goal of SMILE is to ...

## Installation

You can install the development version of SMILE from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JianqiaoWang/SMILE")
```

## Example

The package contains two parts: (1) a general function for SMILE estimation and inference based on the input of (S, y); (2) the function designed for the genetic application. 
The input is a genotype  bed.file, specifed weight matrix W or a block division, and the outcome file . 

# Statistical evaluation

```{r eval = F}
set.seed(123)
library(dplyr)
library(parallel)
n = 800
p = 1000
h2 = 0.4
rho = 0.9
coef_type = "sparse"

# autoregressive structure with rho = 0.9
Sigma = (rho)^(abs(outer(1:p, 1:p, "-")))
Sigma_sqrt = with(eigen(Sigma), vectors %*% (values^(1/2) * t(vectors)))


if(coef_type == "normal"){
beta_coef = rnorm(p)
}


if(coef_type == "sparse"){
  beta_coef = rep( c(1, rep(0, 49)), 20  )
}

X = matrix(rnorm(n*p), n, p) %*% Sigma_sqrt
beta_coef = beta_coef/ c(sqrt(t( beta_coef) %*% Sigma %*% beta_coef )) * sqrt(h2)

# generate the outcome
eps = rnorm(n) * sqrt(1 - h2)
y = X %*% beta_coef + eps

# center y and X and  calculate the heritability

y = scale(y, scale = F)
X = scale(X, scale = F)
scoring_vanila(X, y, W = NULL, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F)
```
Simulate the heterogenous genetic effects 

```{r eval = F}

#ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD <- readRDS("/n/holystore01/LABS/xlin/Lab/wangjq/svd_heri/prepare/ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD.rds")
data('ukbb_geno_EUR_chr1_4.8M_5.5M_LD')
Sigma = ukbb_geno_EUR_chr1_4.8M_5.5M_LD[1:1000, 1:1000]
Sigma_sqrt = with(eigen(Sigma), vectors %*% (values^(1/2) * t(vectors)))


# simulate the weighted normal effects  weight * alpha
weight = sqrt(rowSums(Sigma^2))
alpha  = rnorm(p,  1)
weight[-sample(1:1000, 800)] = 0
beta_coef = weight * alpha
beta_coef = beta_coef/sqrt( as.numeric( t( beta_coef) %*% Sigma %*% beta_coef ) ) * sqrt(h2)


X = matrix(rnorm(n*p), n, p) %*% Sigma_sqrt
# generate the outcome
eps = rnorm(n) * sqrt(1 - h2)
y = X %*% beta_coef + eps


y = scale(y, scale = F)
X = scale(X, scale = F)
scoring_vanila(X, y, W = NULL, theta0 =c(0.5, 0.5), max_iter = 120, eps = 1e-6, verbose = F) # biased estimation with XX\tr


```


# Genetic function and applications 

We only need the bed.file and specified  W matrix
```{r eval = F}
library(bigsnpr)
bedfile = "/n/holystore01/LABS/xlin/Lab/wangjq/UKBB/geno/ukb_cal_autochr_v2_EUR_0_01.bed"
obj.bed <- bed(bedfile)
n = nrow(obj.bed)


## ---------------------------------------------------
weight_mat = readRDS("/n/holystore01/LABS/xlin/Lab/wangjq/UKBB/prepare/weight_mat.rds")
total_rank = sum(unlist(weight_mat[[2]]))
weight_mat = weight_mat[[1]]


## -- match the phenotypes ---------------
pheno = data.table::fread(pheno.file) %>% dplyr::select(FID, IID, all_of(outcome)) %>% dplyr::rename(Y=outcome)
pheno_info = dplyr::left_join(obj.bed$fam, pheno,  by = c( "family.ID" = "FID", "sample.ID" = "IID" ) ) %>% 
  dplyr::mutate(index = 1:n()) ;
pheno_info_non_missing = pheno_info[!is.na(pheno_info$Y), ]


## calculate the scaling
stats_X =bed_scaleBinom(obj.bed, ind.row = pheno_info_non_missing$index, ncores = nb_cores())
means = stats_X$center
sds = stats_X$scale

res = fisher_scoring_var_FBM(obj.bed, scale(pheno_info_non_missing$Y, scale = T), theta0 = theta0, ind.row = pheno_info_non_missing$index,  
                             method = "cg",  weight_mat = weight_mat,  weight_rank = total_rank,  approx = "MC", verbose = T, scale = stats_X$scale*sqrt(total_rank), 
                       nsim = 60, center = stats_X$center, ncores = nb_cores()) # test the bigstatsr pakcgae

```
