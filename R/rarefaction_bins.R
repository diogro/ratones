library(doParallel)
registerDoParallel(cores = 2)

apply(raw.main.data$control.t[,10:44], 1, gm_mean) %>% summary()
main.data$control.t$P49 %>% summary()

fat.t <- raw.main.data$control.t %>% filter(P49 > mean(P49)) %>% makeMainData
light.t <- raw.main.data$control.t %>% filter(P49 < mean(P49)) %>% makeMainData
big.t <- raw.main.data$control.t [apply(raw.main.data$control.t[,10:44], 1, gm_mean) > main.data$control.t$gm_mean,]  %>% makeMainData
small.t <- raw.main.data$control.t [apply(raw.main.data$control.t[,10:44], 1, gm_mean) < main.data$control.t$gm_mean,]  %>% makeMainData

Bins= list("big.t" = big.t,
           "small.t" = small.t,
           "fat.t" = fat.t,
           "light.t" = light.t,
           "control.t" = main.data$control.t,
           "upwards.h'" = main.data$`upwards.h'`,
           "upwards.s'" = main.data$`upwards.s'`,
           "downwards.h" = main.data$downwards.h,
           "downwards.s" = main.data$downwards.s )
bins.MMxStats = Bins %>% llply(., function(x) {CalculateMatrix(lm(as.matrix(dplyr::select(x$full, IS_PM:BA_OPI)) ~ x$full$AGE*x$full$SEX)) } ) %>% ldply(., MeanMatrixStatistics)

raref.R2.Bins = Bins %>% llply(., function(x) {RarefactionStat(x$ed, cor ,function(x, y) CalcR2(y), num.reps = 100, parallel = T, progress = "text")})

raref.r2.t = rbind(raref.R2.Bins$fat.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Heavier.t"),
                       raref.R2.Bins$light.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>%mutate("Bin" = "Lighter.t"),
                       raref.R2.Bins$control.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Control.t"),
                       raref.R2.Bins$`upwards.h'` %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Upwards.h'"),
                       raref.R2.Bins$`upwards.s'` %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Upwards.s'"),
                       raref.R2.Bins$downwards.h %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Downwards.h"),
                       raref.R2.Bins$downwards.s %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Downwards.s") )
raref.r2.t$Bin = factor(raref.r2.t$Bin, levels = unique(raref.r2.t$Bin) )

  raref.r2.t %>% ggplot() +
  geom_linerange(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), ymin = Min., ymax = Max., group = Bin, col = Bin) ) +
  geom_point(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), y= Mean, group = Bin, col = Bin, line = Bin), size = 1.7) +
  scale_colour_manual(values = c("red", "blue", viridis(5))) +
  #facet_wrap(~Bin, ncol = 1) + 
  theme_cowplot() + background_grid(major = 'y', size.major = 0.5, size.minor = 0.4) +
  theme(legend.position = c(0.9, 0.8), panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.4)) + 
  xlab("Sample size") + ylab("Mean squared correlation") + guides(color=guide_legend(title="Line")) + ggtitle("") +
 

save_plot("raref.r2.t.plot.png", raref.r2.t.plot, base_aspect_ratio = 1, base_height = 7, ncol = 1)

raref.PC1.Bins = Bins %>% llply(., function(x) {RarefactionStat(x$ed, cov ,function(x, y) Pc1Percent(y), num.reps = 100, parallel = T)})

raref.pc1.t.Weight = rbind(raref.PC1.Bins$fat.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Heavier.t"),
                       raref.PC1.Bins$light.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>%mutate("Bin" = "Lighter.t"),
                       raref.PC1.Bins$control.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Control.t"))
raref.pc1.t.Weight$Bin = factor(raref.t.Weight$Bin, levels = unique(raref.t.Weight$Bin) )

raref.pc1.t.Size = rbind(raref.PC1.Bins$big.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Larger.t"),
                     raref.PC1.Bins$small.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>%mutate("Bin" = "Smaller.t"),
                     raref.PC1.Bins$control.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Control.t"))
raref.pc1.t.Size$Bin = factor(raref.t.Size$Bin, levels = unique(raref.t.Size$Bin) )

raref.pc1.t.Weight.plot = raref.t.Weight %>% ggplot() +
  geom_linerange(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), ymin = Min., ymax = Max., group = Bin, col = Bin) ) +
  geom_point(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), y= Mean, group = Bin, col = Bin, line = Bin)) +
  scale_colour_manual(values = c("red", "blue", viridis(5)[1])) +
  #facet_wrap(~Bin, ncol = 1) + 
  theme_cowplot() + background_grid(major = 'y', size.major = 0.5, size.minor = 0.4) +
  theme(legend.position = c(0.9, 0.8), panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.4)) + 
  xlab("") + ylab("Mean squared correlation") + guides(color=guide_legend(title="Line")) + ggtitle("Control bins separated by weight") +
  ylim(0.05, 0.2) + xlim(10,60)
raref.pc1.t.Size.plot = raref.t.Size %>% ggplot() +
  geom_linerange(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), ymin = Min., ymax = Max., group = Bin, col = Bin) ) +
  geom_point(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), y= Mean, group = Bin, col = Bin, line = Bin)) +
  scale_colour_manual(values = c("red", "blue", viridis(5)[1])) +
  #facet_wrap(~Bin, ncol = 1) + 
  theme_cowplot() + background_grid(major = 'y', size.major = 0.5, size.minor = 0.4) +
  theme(legend.position = c(0.9, 0.8), panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.4)) + 
  xlab("Sample size") + ylab("") + guides(color=guide_legend(title="Line")) + ggtitle("Control bins separated by cranial size") +
  ylim(0.05, 0.2) + xlim(10, 60)

delta_Z = full_data %>% filter(selection == "upwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans - full_data %>% filter(selection == "downwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans
Z_control <- full_data %>% filter(selection == "control") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans 
delta_Z / Z_control 
