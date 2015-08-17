source("R/run_ratones_MCMCglmm_P.R")

vectorCor <- function(x, y) t(Normalize(x)) %*% Normalize(y)

ML_R2 = ldply(main.data, function(x) CalcR2(x$cov.matrix))

r_models %>% llply(function(x) x$Ps) %>% ldply(function(x) adply(x, 1, CalcR2)) %>%
  group_by(.id) %>% summarise_each(., funs(mean, find_CI_lower, find_CI_upper), V1) %>%
  ggplot(aes(.id, mean, group = .id)) + geom_point() + geom_errorbar(aes(ymin = find_CI_lower, 
                                                                         ymax = find_CI_upper)) +
  geom_point(data = ML_R2, color = 'red', aes(.id, V1)) + 
  theme_bw() + labs(y = expression(R^2))

#mcmc_stats = tbl_df(ldply(r_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics), .parallel = TRUE))
#save(mcmc_stats, file = "./Rdatas/mcmc_stats")
load("./Rdatas/mcmc_stats")
names(mcmc_stats)[4] <- 'pc1.percent'
mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>%
  separate(.id, c( 'treatment', 'strain')) %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() +
  facet_wrap(~variable, scale = 'free') + theme_bw() + scale_fill_manual(values = c(c, h, s)) -> global_stats_plot

x = main.data[[2]]

traits = select(full_data, P49, IS_PM:BA_OPI)
traits$P49 %<>% log
delta_Zs <- llply(main.data, function(x) x$ed.means - main.data$control.control$ed.means)
delta_Zs <- llply(main.data, function(x) x$ed.means - colMeans(traits))
plsr <- llply(main.data, function(x) x$plsr)

directionalVariation <- function(cov.matrix, strain){
  delta_Z <- delta_Zs[[strain]]
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(cov.matrix, ret.dim = 13)[[1]], delta_Z)  
  beta = c(1, rep(0, 35))
  data.frame(corDZDZ = abs(vectorCor(cov.matrix %*% beta, delta_Z)),
             DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])),
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

DzPC1 <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'DZpc1') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + theme_bw() + 
  ggtitle("Correlation of mean change and PC1") + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation")

corDZDZ <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'corDZDZ') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0.32) +
  ggtitle("Correlation of observed change\n and expected change in multivariate mean") + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation")

treatment %<>%  separate(.id, c('treatment', 'strain'))
normDZ_DzPC1 = ggplot(treatment, aes(normDZ, DZpc1, group = interaction(treatment, strain), color = strain)) + geom_violin(aes(fill = strain), alpha = 0.3) + geom_jitter(aes(shape = treatment), size = 3, position = position_jitter(width = .03)) + scale_fill_manual(values = c(h, s)) + scale_color_manual(values = c(h, s))

