if(!require (mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if(!require (ggplot2)) {install.packages("ggplot2"); library(ggplot2)}

source("R/read_ratones.R")

beta = rep(1, 35)

simulateBetaEstimation = function(beta_center, G, beta_var = 0.5) {
  beta = beta_center + rnorm(35, 0, beta_var)
  delta_Z = G %*% beta
  pop = mvrnorm(60, mu = rep(0, 35), Sigma = G)
  G_estimated = cov(pop)
  beta_estimated = solve(G_estimated, delta_Z)
  beta_non_noise = solve(ExtendMatrix(G_estimated, ret.dim = 10)[[1]], delta_Z)
  data.frame(beta = vectorCor(beta_estimated, beta),
             beta_non_noise = vectorCor(beta_non_noise, beta),
             delta_z = vectorCor(delta_Z, beta))
}
result = rdply(1000, simulateBetaEstimation(beta, delta_Z)) 
m_result = gather(result, estimator, value, beta:delta_z)
figure_S7 = ggplot(m_result, aes(value, fill = estimator)) + 
  geom_density() + labs(x= expression(paste("Vector correlation with the real ", beta)))+
  scale_fill_brewer(palette = "Greys", "Estimators", labels = c(expression(paste(beta, " - from equation" ), paste( beta, " - non noise"), paste(Delta, "z")) ),
  guide = guide_legend(
    label.position = "right",
    label.hjust = 0))
figure_S7  
save_plot("figureS7.png", figure_S7, base_aspect_ratio = 1.3, base_height = 4)
