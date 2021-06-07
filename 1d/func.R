library(MASS)

b <- function(y, t, alpha, sigma, theta){
  b_de <- rep(0, 2)
  b_nu <- rep(0, 2)
  for (i in 1:2){
    alp <- alpha[i]
    sig <- sigma[i]
    
    sig_t <- (1-t+t*sig)
    b_de[i] <- exp(-(t*alp^2+y^2-y^2*sig-2*y*alp)/(2*sig_t))/(sig_t^0.5)
    b_nu[i] <- b_de[i]*(alp+y*(sig-1))/sig_t
  }
  return (b_nu %*% theta/theta%*%b_de)
}

sample_1d<-function(K, alpha, sigma, theta){
  y <- matrix(0, K)
  s <- 1/K
  for (k in 1:(K-1)){
    eps <- rnorm(1, 0, 1)
    y[(k+1)] <- y[k] + s*b(y[k], k*s, alpha, sigma, theta) +  sqrt(s)*eps
  }
  return(y[K])
}


SFS_1d<-function(N, K, alpha, sigma, theta, parallel=FALSE){
  if (parallel==FALSE){
    res <- matrix(0, N)
    for (i in 1:N){
      res[i] <- sample_1d(K, alpha, sigma, theta)
    }
  }
  else {
    library(doParallel)
    cl <- makeCluster(parallel)
    registerDoParallel(cl)
    res <- foreach(x=1:N, .combine='rbind', .packages='MASS', .export=c("sample_1d", "b")) %dopar% sample_1d(K, alpha, sigma, theta)
    stopCluster(cl)
  }
  return(res)
}
