source("R/run_ratones_MCMCglmm_P.R")
ML_R2 = ldply(main.data, function(x) CalcR2(x$cov.matrix))

r_models %>% llply(function(x) x$Ps) %>% ldply(function(x) adply(x, 1, CalcR2)) %>%
  group_by(.id) %>% summarise_each(., funs(mean, find_CI_lower, find_CI_upper), V1) %>%
  ggplot(aes(.id, mean, group = .id)) + geom_point() + geom_errorbar(aes(ymin = find_CI_lower, 
                                                                         ymax = find_CI_upper)) +
  geom_point(data = ML_R2, color = 'red', aes(.id, V1)) + 
  theme_bw() + labs(y = expression(R^2))

mcmc_stats = tbl_df(ldply(r_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics), .parallel = TRUE))
mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>%
  separate(.id, c( 'treatment', 'strain')) %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, variable), color = strain)) + geom_boxplot() +
  facet_wrap(~variable, scale = 'free') + theme_bw()

x = main.data[[2]]
delta_Zs <- llply(main.data, function(x) x$ed.means - main.data$control.control$ed.means)
plsr <- llply(main.data, function(x) x$plsr)

directionalVariation <- function(cov.matrix, strain){
  delta_Z <- delta_Zs[[strain]]
  plsr_strain <- plsr[[strain]]
  plsr_control <- plsr[[1]]
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(cov.matrix, ret.dim = 13)[[1]], delta_Z)  
  beta = c(1, rep(0, 35))
  #eval1 <- eigen(cov.matrix)$value[1]
  #eval1_ <- (eigen(cov.matrix[-1, -1])$value)
  #mat.stat <- MeanMatrixStatistics(cov.matrix, parallel = TRUE)
  data.frame(corDZDZ = vectorCor(cov.matrix %*% beta, delta_Z),
             flex_DZ = Flexibility(cov.matrix, Normalize(delta_Z)),#/mat.stat['flexibility'],
             flex_beta = Flexibility(cov.matrix, Normalize(beta)),#/mat.stat['flexibility'], 
             flex_plsr_strain = Flexibility(cov.matrix[-1,-1], Normalize(plsr_strain)),#/mat.stat['flexibility'], 
             flex_plsr_control = Flexibility(cov.matrix[-1,-1], Normalize(plsr_control)),#/mat.stat['flexibility'], 
             Evol_DZ = Evolvability(cov.matrix, Normalize(delta_Z)),#/eval1,
             #DZpc1 = abs(vectorCor(delta_Z,
             #                       eigen(cov.matrix)$vector[,1])),
             #corBetaBetaS = vectorCor(beta, beta_s),
             #corBetaBetaNN = vectorCor(beta, beta_NN),
             #corDZbetaS = vectorCor(delta_Z, beta_s),            
             Evol_beta = Evolvability(cov.matrix, Normalize(beta)),#/eval1, 
             
             Evol_plsr_strain = Evolvability(cov.matrix[-1, -1], Normalize(plsr_strain)),#/eval1_, 
             
             Evol_plsr_control = Evolvability(cov.matrix[-1, -1], Normalize(plsr_control)),#/eval1_, 
             normDZ = Norm(delta_Z))
}

treatment <- ldply(r_models[-1], function(model) adply(model$Ps, 1, 
                                                       directionalVariation, 
                                                       model$strain), .parallel = TRUE)
treatment$type <- 'treatment'
control   <- ldply(r_models[-1], function(model) adply(r_models[['control.control']]$Ps, 1, 
                                                       directionalVariation, 
                                                       model$strain), .parallel = TRUE)
control$type <- 'control'
stats <- melt(rbind(treatment, control))[-2]
stats$fe = '0'
stats$fe[grep('flex', as.character(stats$variable))] <- 'flex'
stats$fe[grep('Evol', as.character(stats$variable))] <- 'Evol'

stats %>% filter(variable != 'normDZ') %>% 
  ggplot(., aes(variable, value, group = interaction(type, variable, .id), color = type)) + geom_boxplot() + facet_grid(.id~fe, scales = "free")

stats %>% filter(variable != 'normDZ') %>% 
  ggplot(., aes(.id, value, group = interaction(type, variable, .id), color = type)) + geom_boxplot() + facet_wrap(~variable, scales = "free")
