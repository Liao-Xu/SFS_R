b <- function(y, p, kap, t, alpha, sigma, theta, stats){
  Ip <- diag(p)
  b_de <- rep(0, kap)
  b_nu <- matrix(0, p, kap)
  for (i in 1:kap){
    alp <- t(alpha[i,, drop = FALSE])
    sig <- sigma[,,i]
    
    inv_sig <- stats[[1]][[i]]                     # sigma^(-1)
    is_a <- stats[[2]][[i]]                        # sigma^(-1)*alpha
    at_is_a <- stats[[3]][[i]]                     # alpha^T*sigma^(-1)*alpha
    t_is <- (1-t)*inv_sig                          # (1-t)*sigma
    t_is_y <- (t_is %*% alp + y)                   # (1-t)*sigma*alpha + y
    com_I <-solve(t*Ip + t_is)                     # ((1-t)*sigma+t*I)^(-1)
    com_sig <- det(t*sig + (1-t)*Ip)^0.5           # (t*sigma+(1-t)*I)^(-0.5)
    
    g <- exp((sum((com_I^0.5 %*% (t_is%*%alp+y))^2)-sum(y^2))/(2-2*t) - 0.5*at_is_a)
    b_de[i] <- g/com_sig
    b_nu[,i] <- (is_a + (Ip-inv_sig) %*% com_I %*% t_is_y)*b_de[i]
  }
  return (drop(b_nu %*% theta)/sum(theta*b_de))
}

sample_2d <- function(p, kap, K, alpha, sigma, theta, stats){
  y <- matrix(0, p, K)
  s <- 1/K
  for (k in 1:(K-1)){
    eps <- mvrnorm(1, rep(0,p), diag(p))
    y[,(k+1)] <- y[,k] + s*b(y[,k], p, kap, k*s, alpha, sigma, theta, stats) +  sqrt(s)*eps
  }
  return(y[,K])
}

SFS_2d<-function(N, K, p, kap, alpha, sigma, theta, parallel=FALSE){
  stats1 <- vector(mode = "list", length = kap)
  stats2 <- vector(mode = "list", length = kap)
  stats3 <- vector(mode = "list", length = kap)
  for (i in 1:kap){
    alp <- t(alpha[i,, drop = FALSE])
    sig <- sigma[,,i]
    
    stats1[[i]] <- solve(sig)                          # sigma^(-1)
    stats2[[i]] <- stats1[[i]] %*% alp                 # sigma^(-1)*alpha
    stats3[[i]] <- t(alp) %*% stats2[[i]]              # alpha^T*sigma^(-1)*alpha
  }
  stats <- list(stats1, stats2, stats3)
  
  if (parallel==FALSE){
    res <- matrix(0, p, N)
    for (i in 1:N){
      res[, i] <- sample_2d(p, kap, K, alpha, sigma, theta, stats)
    }
  }
  else {
    library(doParallel)
    cl <- makeCluster(parallel)
    registerDoParallel(cl)
    res <- foreach(x=1:N, .combine='rbind', .packages='MASS', .export=c("sample_2d", "b")) %dopar% sample_2d(p, kap, K, alpha, sigma, theta, stats)
    res <- t(res)
    stopCluster(cl)
  }
  return(res)
}