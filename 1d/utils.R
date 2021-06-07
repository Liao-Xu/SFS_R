p_density <- function(x, alpha, sigma, theta){
  p1 = 1/(sqrt(2*pi*sigma[1]))*exp(-(x-alpha[1])^2/(2*sigma[1]))
  p2 = 1/(sqrt(2*pi*sigma[2]))*exp(-(x-alpha[2])^2/(2*sigma[2]))
  p3 = theta[1]*p1+(1-theta[1])*p2
  return(p3)
}

plots_1d <- function(res, N, K, alpha, sigma, theta, resDir = "results", bw = "SJ"){
  mainDir = getwd()
  resPath = file.path(mainDir, resDir)
  
  # Create dir for saving results
  ifelse(!dir.exists(resPath), dir.create(resPath), FALSE)
  
  save(res, file=file.path(resPath, paste0("1d_K_",K,"_N_",N,".rdata")))
  x0=seq(alpha[1]-2,alpha[2]+2,by=0.01)
  xnorm=p_density(x0, alpha, sigma, theta)
  
  d <- density(res, bw = bw)
  png(file=file.path(resPath,  paste0("1d_K_", K,"_N_",N,"_bw_",bw,".png")), width = 900, height = 600, res=150)
  plot(d, ylim=c(0,max(xnorm,d$y)), col="red" , lwd=0.8, main="")
  polygon(x0, xnorm,col=gray(0.5), border=0)
  lines(d, col="red", lwd=0.8)
  dev.off()
}