### This repository is the R demo implementation of Schrodinger-Follmer Sampler for Gaussian Mixtures.
### Input: N, K, alpha (kap vector), sigma (kap  vector), theta (kap vector)

# Install missing packages
rm(list=ls(all=TRUE));
list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load local functions
source('func.R')
source('utils.R')


start.time <- Sys.time()
### ======================= initializing ==================
N <- 5000
K <- 200
kap <- 2
alpha <- c(-4, 4)
sigma <- c(0.03, 0.03)
theta <- c(0.5, 0.5)

### ======================= sampling ==================
res <- SFS_1d(N, K, alpha, sigma, theta, parallel=6)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

plots_1d(res, N, K, alpha, sigma, theta, resDir = "results")
