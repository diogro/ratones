source("R/run_ratones_MCMCglmm_P.R")
ML_R2 = ldply(main.data, function(x) CalcR2(x$cov.matrix))
r_models %>% llply(function(x) x$Ps) %>% ldply(function(x) adply(x, 1, CalcR2)) %>%
  group_by(.id) %>% summarise_each(., funs(mean, find_CI_lower, find_CI_upper), V1) %>%
  ggplot(aes(.id, mean, group = .id)) + geom_point() + geom_errorbar(aes(ymin = find_CI_lower, 
                                                                         ymax = find_CI_upper)) +
  geom_point(data = ML_R2, color = 'red', aes(.id, V1)) + 
  theme_bw() + labs(y = expression(R^2))

#mcmc_stats = tbl_df(ldply(r_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics), .parallel = TRUE))
mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>%
  separate(.id, c( 'strain', 'treatment')) %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, variable), color = strain)) + geom_boxplot() +
  facet_wrap(~variable, scale = 'free') + theme_bw()

