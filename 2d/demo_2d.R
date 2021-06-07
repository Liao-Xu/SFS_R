### This repository is the R demo implementation of Schrodinger-Follmer Sampler for Gaussian Mixtures.
### Input: N, K, alpha (p x kap matrix), sigma (p x p x kap  matrix), theta (kap vector)

# Install missing packages
rm(list=ls(all=TRUE));
list.of.packages <- c("MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load local functions
source('func.R')
source('utils.R')

start.time <- Sys.time()
### ======================= initializing ==================
N <- 5000
K <- 100
kap <- 16
rad <- 8 
p <- 2
alpha=matrix(0, kap, p)
for(i in 1:kap){
  alpha[i,1] <- rad * cos(i*2*pi/(kap))
  alpha[i,2] <- rad * sin(i*2*pi/(kap))
}
sigma <- replicate(kap, 0.03*diag(p), simplify="array")
theta <- rep(1/kap, kap)

### ======================= sampling ==================
res <- SFS_2d(N, K, p, kap, alpha, sigma, theta, parallel=6)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
plots_2d(res, N, K, alpha, sigma, theta, resDir = "results")
