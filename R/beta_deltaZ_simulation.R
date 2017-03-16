library(mvtnorm)
library(ggplot2)

source("R/read_ratones.R")

beta = rep(1, 35)

simulateBetaEstimation = function(beta_center, G, beta_var = 0.5) {
  beta = beta_center + rnorm(35, 0, beta_var)
  delta_Z = G %*% beta
  pop = mvrnorm(60, mu = rep(0, 35), Sigma = G)
  G_estimated = cov(pop)
  beta_estimated = solve(G_estimated, delta_Z)
  #beta_non_noise = solve(ExtendMatrix(G_estimated)[[1]], delta_Z)
  data.frame("estimated beta" = vectorCor(beta_estimated, beta),
   #         beta_non_noise = vectorCor(beta_non_noise, beta),
             delta_z = vectorCor(delta_Z, beta))
}
result = rdply(1000, simulateBetaEstimation(beta, P, 0.8)) 
m_result = gather(result, estimator, value, estimated.beta:delta_z)
ggplot(m_result, aes(value, fill = estimator)) + geom_density(alpha = 0.7) + labs(x = "Correlation with true selection gradient", y = "Density")
