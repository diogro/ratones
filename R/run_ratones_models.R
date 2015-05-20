source('./R/read_ratones.R')

vectorCor <- function(x, y) t(Normalize(x)) %*% Normalize(y)

stats <- llply(main.data, function(x) x$cov.matrix) %>% ldply(MeanMatrixStatistics)
ggplot(stats, aes(MeanSquaredCorrelation, flexibility)) + geom_text(aes(label = .id)) + geom_smooth(method = 'lm')

ddply(raw.data, .(strain, treatment), function(x) sd(x$P49)/mean(x$P49))

r2.df <- llply(main.data, function(x) x$cov.matrix) %>% ldply(MonteCarloR2, 50) %>% gather(.id)
names(r2.df)[1] <- 'strain'
r2.df %>% separate(strain, c( 'treatment', 'strain')) %>%
  ggplot(., aes(treatment, value, group = interaction(treatment, strain), color = strain)) +
  geom_boxplot() + 
  theme_classic() + labs(y = expression(R^2))

MatrixCompare(cov2cor(main.data[[1]]$cov.matrix), 
              cov2cor(main.data[[3]]$cov.matrix))

delta_beta_cor <- function(x){
  delta_Z <- x$ed.means - main.data[[1]]$ed.means
  beta_s <- solve(x$cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(x$cov.matrix, ret.dim = 13)[[1]], delta_Z)  
  beta = c(1, rep(0, 35))
  data.frame(corDZDZ = vectorCor(x$cov.matrix %*% beta, delta_Z),
             flex_DZ = Flexibility(Normalize(delta_Z), x$cov.matrix),
             Evol_DZ = Evolvability(Normalize(delta_Z), x$cov.matrix),
             DZpc1_I = abs(vectorCor(eigen(main.data[[1]]$cov.matrix)$vectors[,1],
                                delta_Z)),
             DZpc1 = abs(vectorCor(delta_Z,
                                eigen(x$cov.matrix)$vector[,1])),
             #corBetaBetaS = vectorCor(beta, beta_s),
             #corBetaBetaNN = vectorCor(beta, beta_NN),
             #corDZbetaS = vectorCor(delta_Z, beta_s),
             flex_beta = Flexibility(Normalize(beta), x$cov.matrix), 
             Evol_beta = Evolvability(Normalize(beta), x$cov.matrix), 
             pc1_DZ = Evolvability(Normalize(delta_Z), x$cov.matrix),
             normDZ = Norm(delta_Z))
}
ldply(main.data, delta_beta_cor)
select(stats, .id, flexibility, evolvability)


