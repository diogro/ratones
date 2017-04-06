
evolBeta <- stats %>% dplyr::select(-X1) %>% ggplot(aes(evolBeta, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Scaled directional evolvability (Beta)", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))
condevolBeta <- stats %>% dplyr::select(-X1) %>% ggplot(aes(condevolBeta, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Scaled directional condititonal evolvability (Beta)", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))
BetaPC1 <- stats %>% dplyr::select(-X1) %>% ggplot(aes(BetaPC1, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = expression(paste("Vector correlation of ", beta, " and the first eigenvector")), y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(text = element_text(size = 20))

figure_4_Beta <- ggdraw() +
  draw_plot(evolBeta, 0, 0.5, 0.5, 0.5) +
  draw_plot(condevolBeta, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(BetaPC1, 0, 0, 0.7, 0.5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0), c(1, 1, 0.5), size = 20)
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure4_beta.png", figure_4_Beta, ncol = 2, nrow = 2, base_aspect_ratio = 1.3, base_height = 4)


betaDzComparison <- function(cov.matrix, line){
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(cov.matrix, ret.dim = 13)[[1]], delta_Z)
  data.frame(beta  = vectorCor(delta_Z, beta_s), 
             betaN = vectorCor(delta_Z, beta_NN))
}

stats_beta <- ldply(g_models[2:5], function(model) adply(model$Ps, 1,
                                                         betaDzComparison,
                                                         model$line), .parallel = TRUE)

DzBeta <- stats_beta %>% dplyr::select(-X1) %>% ggplot(aes(beta, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Beta x Delta Z", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(text = element_text(size = 20))
DzBetaN <- stats_beta %>% dplyr::select(-X1) %>% ggplot(aes(betaN, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Beta x Delta Z", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))
