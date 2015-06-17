source("R/run_ratones_MCMCglmm_P.R")

vectorCor <- function(x, y) t(Normalize(x)) %*% Normalize(y)

ML_R2 = ldply(main.data, function(x) CalcR2(x$cov.matrix))

r_models %>% llply(function(x) x$Ps) %>% ldply(function(x) adply(x, 1, CalcR2)) %>%
  group_by(.id) %>% summarise_each(., funs(mean, find_CI_lower, find_CI_upper), V1) %>%
  ggplot(aes(.id, mean, group = .id)) + geom_point() + geom_errorbar(aes(ymin = find_CI_lower, 
                                                                         ymax = find_CI_upper)) +
  geom_point(data = ML_R2, color = 'red', aes(.id, V1)) + 
  theme_bw() + labs(y = expression(R^2))

# mcmc_stats = tbl_df(ldply(r_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics), .parallel = TRUE))
# save(mcmc_stats, file = "./Rdatas/mcmc_stats")
load("./Rdatas/mcmc_stats")
mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>%
  separate(.id, c( 'treatment', 'strain')) %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() +
  facet_wrap(~variable, scale = 'free') + theme_bw() -> global_stats_plot

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
  flexibility <- Flexibility(cov.matrix)
  evolvability <- Evolvability(cov.matrix)
  data.frame(corDZDZ = vectorCor(cov.matrix %*% beta, delta_Z),
             corDZplsrDZ = vectorCor(cov.matrix[-1,-1] %*% plsr_strain, delta_Z[-1]),
             flex_DZ = Flexibility(cov.matrix, Normalize(delta_Z))/flexibility,
             flex_beta = Flexibility(cov.matrix, Normalize(beta))/flexibility, 
             #flex_plsr_strain = Flexibility(cov.matrix[-1,-1], Normalize(plsr_strain)),#/mat.stat['flexibility'], 
             #flex_plsr_control = Flexibility(cov.matrix[-1,-1], Normalize(plsr_control)),#/mat.stat['flexibility'], 
             Evol_DZ = Evolvability(cov.matrix, Normalize(delta_Z))/evolvability,
             DZpc1 = abs(vectorCor(delta_Z,
                                    eigen(cov.matrix)$vector[,1])),
             #corBetaBetaS = vectorCor(beta, beta_s),
             #corBetaBetaNN = vectorCor(beta, beta_NN),
             #corDZbetaS = vectorCor(delta_Z, beta_s),            
             Evol_beta = Evolvability(cov.matrix, Normalize(beta))/evolvability
             
             #Evol_plsr_strain = Evolvability(cov.matrix[-1, -1], Normalize(plsr_strain)),#/eval1_, 
             
             #Evol_plsr_control = Evolvability(cov.matrix[-1, -1], Normalize(plsr_control)),#/eval1_, 
             #normDZ = Norm(delta_Z))
  )
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

DzPC1 <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'DZpc1') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + theme_bw() + ggtitle("Correlation of mean change and PC1")

testplsrHip <- function (x) {
  hips = matrix(c(as.numeric(x$ed.means - main.data[[1]]$ed.means > 0), 
                  as.numeric(x$ed.means - main.data[[1]]$ed.means < 0)), ncol = 2)
  TestModularity(cov2cor(x$cov.matrix), hips)
}
ldply(main.data, testplsrHip)
LModularity((main.data[[1]]$cov.matrix))
