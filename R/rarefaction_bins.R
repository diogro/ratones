library(doParallel)
registerDoParallel(cores = 2)

fat.t <- raw.main.data$control.t %>% filter(P49 > mean(P49)) %>% makeMainData
light.t <- raw.main.data$control.t %>% filter(P49 < mean(P49)) %>% makeMainData

unique(fat.t$info$ID) %>% length()
unique(light.t$info$ID) %>% length()

CalculateMatrix(lm(as.matrix(dplyr::select(fat.t$full, IS_PM:BA_OPI)) ~ fat.t$full$AGE*fat.t$full$SEX)) %>% MeanMatrixStatistics()
CalculateMatrix(lm(as.matrix(dplyr::select(light.t$full, IS_PM:BA_OPI)) ~ light.t$full$AGE*light.t$full$SEX)) %>% MeanMatrixStatistics()
CalculateMatrix(lm(as.matrix(dplyr::select(main.data$control.t$full, IS_PM:BA_OPI)) ~ main.data$control.t$full$AGE*main.data$control.t$full$SEX)) %>% MeanMatrixStatistics()

library(doParallel)
registerDoParallel(cores = 2)

raref.fat.t <- RarefactionStat(fat.t$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.light.t <- RarefactionStat(light.t$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.all.t <- RarefactionStat(main.data$control.t$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.all.s <- RarefactionStat(main.data$downwards.s$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.all.s.up <- RarefactionStat(main.data$`upwards.s'`$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.all.h <- RarefactionStat(main.data$downwards.h$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)
raref.all.h.up <- RarefactionStat(main.data$`upwards.h'`$ed, cor, function(x, y) CalcR2(y), num.reps = 100, parallel = T)



raref.t = rbind(raref.fat.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Heavier.t"),
raref.light.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>%mutate("Bin" = "Lighter.t"),
raref.all.t %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Control.t"), 
raref.all.h.up %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Upwards.h'"),
raref.all.s.up %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Upwards.s'"),
raref.all.h %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Downwards.h"),
raref.all.s %>% ldply(summary) %>% filter(as.numeric_version(X1) > 10) %>% mutate("Bin" = "Downwards.s") )

raref.t$Bin = factor(raref.t$Bin, levels = unique(raref.t$Bin) )

raref.t %>% ggplot() +
  geom_linerange(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), ymin = Min., ymax = Max., group = Bin, col = Bin) ) +
  geom_point(aes(x = (as.numeric(X1) + as.numeric(Bin) *0.1), y= Mean, group = Bin, col = Bin, line = Bin)) +
  scale_colour_manual(values = c("red", "blue", viridis(5))) +
  #facet_wrap(~Bin, ncol = 1) + 
  theme_cowplot() + background_grid(major = 'y', size.major = 0.5, size.minor = 0.4) +
  theme(legend.position = c(0.8, 0.8), panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.4)) + 
  xlab("Sample size") + ylab("Mean squared correlation") + guides(color=guide_legend(title="Line"))


delta_Z = full_data %>% filter(selection == "upwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans - full_data %>% filter(selection == "downwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans
Z_control <- full_data %>% filter(selection == "control") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans 
delta_Z / Z_control 
