## code to prepare `LD` dataset goes here
ukbb_geno_EUR_chr1_4.8M_5.5M_LD <- readRDS("C:/Users/97683/Dropbox/Git-rep/SMILE/data-raw/ukbb_geno_subset_EUR_chr1_4.8M_5.5M_LD.rds")
usethis::use_data(ukbb_geno_EUR_chr1_4.8M_5.5M_LD, overwrite = TRUE)
