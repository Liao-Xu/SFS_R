library(tidyverse)

plots_2d <- function(res, N, K, alpha, sigma, theta, resDir = "results"){
  mainDir = getwd()
  resPath = file.path(mainDir, resDir)
  
  # Create dir for saving results
  ifelse(!dir.exists(resPath), dir.create(resPath), FALSE)
  
  data <- data.frame(x=res[1,], y=res[2,])
  save(data, file=file.path(resPath, paste0("2d_K_",K,"_N_",N,".rdata")))
  png(file=file.path(resPath,   paste0("2d_K_", K,"_N_",N,"_raster.png")), width = 900, height = 800, res=150)
  p <- ggplot(data, aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_continuous(type = "viridis")+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) 
  print(p)
  dev.off()
  
  png(file=file.path(resPath,  paste0("2d_K_", K,"_N_",N,"_scatter.png")), width = 900, height = 800, res=150)
  p <- ggplot(data, aes(x=x, y=y) ) +
    geom_density_2d()+
    geom_point(alpha=0.2)
  print(p)
  dev.off()
}