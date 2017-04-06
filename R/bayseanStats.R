if(!require(doMC)) {install.packages('doMC'); library(doMC)}
if(!require(viridis)) install.packages("viridis")
library(viridis)
registerDoMC(3)

#source("R/run_ratones_MCMCglmm_P.R")
source("R/read_ratones.R")
p_mcmc_stats = tbl_df(ldply(g_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE))
g_mcmc_stats = tbl_df(ldply(g_models, function(x) adply(x$Gs, 1, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE))g_mcmc_stats = tbl_df(ldply(g_models, function(x) adply(x$Gs, 1, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE))

mcmc_stats = g_mcmc_stats
# save(mcmc_stats, file = "./Rdatas/P_mcmc_stats")
load("./Rdatas/P_mcmc_stats")
names(mcmc_stats) <- gsub("pc1%", "pc1.percent", names(mcmc_stats))

global_stats <- mcmc_stats %>% dplyr::select(.id, MeanSquaredCorrelation, flexibility, pc1.percent, evolvability, conditional.evolvability) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
global_stats$line <- gsub('control', 't', global_stats$line)
global_stats$line <- factor(global_stats$line, levels = lines)

r2_plot <- ggplot(filter(global_stats, variable == "MeanSquaredCorrelation"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border() + labs(x = "Mean squared correlation", y = "Density") + theme(legend.position = c(0.7, 0.7), text = element_text(size = 20))
pc1.percent_plot <- ggplot(filter(global_stats, variable == "pc1.percent"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line") + background_grid(major = 'x', minor = "none") +  panel_border() + labs(x = "Proportion of variation in E1", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
flexibility_plot <- ggplot(filter(global_stats, variable == "flexibility"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line") + background_grid(major = 'x', minor = "none") +  panel_border() + labs(x = "Mean flexibility", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
evolvability_plot <- ggplot(filter(global_stats, variable == "evolvability"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line")+ background_grid(major = 'x', minor = "none") +  panel_border() + labs(x = "Mean evolvability", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
cond_evolvability_plot <- ggplot(filter(global_stats, variable == "conditional.evolvability"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line")+ background_grid(major = 'x', minor = "none") +  panel_border() + labs(x = "Mean conditional evolvability", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))

figure_3 <- ggdraw() +
  draw_plot(r2_plot, 0, 0.5, 0.5, 0.5) +
  draw_plot(pc1.percent_plot, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(flexibility_plot, 0, 0, 0.5, 0.5) +
  draw_plot(evolvability_plot, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 20)
save_plot("figure3.png", figure_3, ncol = 2, nrow = 2, base_aspect_ratio = 1.3, base_height = 4)

directionalVariation <- function(cov.matrix, line){
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(P)[[1]], delta_Z)
  condEvol = (sum(ConditionalEvolvability(cov.matrix)))
  data.frame(evolDZ = Evolvability(cov.matrix, Normalize(delta_Z)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             evolBeta = Evolvability(cov.matrix, Normalize(beta_NN)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             condevolDZ = ConditionalEvolvability(cov.matrix, Normalize(delta_Z)) / condEvol,
             condevolBeta = ConditionalEvolvability(cov.matrix, Normalize(beta_NN)) / condEvol,
             BetaPC1 = abs(vectorCor(beta_NN, eigen(cov.matrix)$vector[,1])),
             DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])))
}

p_stats <- ldply(g_models, function(model) adply(model$Ps, 1,
                                               directionalVariation,
                                               model$line), .parallel = TRUE)
g_stats <- ldply(g_models, function(model) adply(model$Gs, 1,
                                                 directionalVariation,
                                                 model$line), .parallel = TRUE)
stast = g_stats

DzPC1 <- stats %>% dplyr::select(-X1) %>% ggplot(aes(DZpc1, group = .id, fill = .id)) + geom_density(alpha = 0.5) + labs(x = expression(paste("Vector correlation of ", Delta, "z and E1")), y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = c(1.25, 0.5), text = element_text(size = 20))
evolDZ <- stats %>% dplyr::select(-X1) %>% ggplot(aes(evolDZ, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Scaled directional\n evolvability", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))
condevolDZ <- stats %>% dplyr::select(-X1) %>% ggplot(aes(condevolDZ, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Scaled directional\n conditional evolvability", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'x', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))

figure_4 <- ggdraw() +
  draw_plot(evolDZ, 0, 0.5, 0.5, 0.5) +
  draw_plot(condevolDZ, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(DzPC1, 0.2, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0.2), c(1, 1, 0.5), size = 18)

save_plot("figure4_SIversion.png", figure_4, ncol = 2, nrow = 2, base_aspect_ratio = 1.3, base_height = 4)

PCones <- t(laply(g_models, function(x) eigen(x$P)$vectors[,1]))
colnames(PCones) <- c("Control t", "Downwards h", "Downwards s", "Upwards h'", "Upwards s'")
rownames(PCones) <- colnames(main.data[[1]]$ed)

library(xtable)
print.xtable(xtable(PCones, digits = 3)) 

PC1_iso_cor = aaply(PCones, 2, function(x) abs(vectorCor(x, sqrt(rep(1/35, 35)))))
PC1_dz_cor = aaply(PCones, 2, function(x) abs(vectorCor(x, delta_Z)))
vectorCor(sqrt(rep(1/35, 35)), delta_Z)

G_comp = cbind(RandomSkewers(llply(g_models, `[[`, "P"), G)[,1:2],
                      KrzCor(llply(g_models, `[[`, "P"), G)[,2])
G_comp[6,] = data.frame("Within Groups P-matrix", RandomSkewers(P, G)[1], KrzCor(P, G)[1])
G_comp[6,1] = "Within Groups P-matrix"
names(G_comp) = c(".id", "Random Skewers", "Krzanowski")
rownames(G_comp) = NULL
xtable(G_comp)

GP = RandomSkewers(Gs, Ps, parallel = TRUE)

ic.min<-function(x) {
  t<- length(x)
  ic<- sort(x)[round(t*c(0.025) )]
  names(ic) <- c("min")
  return(ic)
}

ic.max<-function(x) {
  t<- length(x)
  ic<-sort(x)[round(t*c(0.975) )] 
  names(ic) <- c("max")
  return(ic)
}

global_stats %>% group_by(variable, selection, line) %>% summarise_each(funs(ic.min, mean, ic.max), value) %>% xtable %>% print.xtable(type = "latex", include.rownames = FALSE)
stats %>% dplyr::select(.id, evolDZ, condevolDZ, DZpc1) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.") %>% group_by(variable, selection, line) %>% summarise_each(funs(ic.min, mean, ic.max), value) %>% xtable() %>% print.xtable(type = "latex", include.rownames = FALSE)

m = g_models$control.t$Ps

ic.mean.mx.dist = function(m) 
  {
CI.min = matrix(NA, nrow = dim(m)[2], ncol = dim(m)[3]) 
CI.max = matrix(NA, nrow = dim(m)[2], ncol = dim(m)[3]) 
Mean = matrix(NA, nrow = dim(m)[2], ncol = dim(m)[3]) 
for (j in 1:dim(m)[3]) {
  for (i in 1:dim(m)[2]) { 
    CI.min[i,j] = ic.min (m[ , i, j])
    CI.max[i,j] = ic.max (m[ , i, j])
    Mean[i,j] = mean(m[ , i, j])
    }
}
return(list("CI.min" = CI.min,
            "CI.max" = CI.max,
            "Mean" = Mean ))
}

ic.mean.P.matrices = g_models %>% llply(., function(x) x$Ps ) %>% llply(., ic.mean.mx.dist)
ic.mean.G.matrices = g_models %>% llply(., function(x) x$Gs ) %>% llply(., ic.mean.mx.dist)

all.lines.P.matrices = g_models %>% llply(., function(x) x$Ps )
all.lines.G.matrices = g_models %>% llply(., function(x) x$Gs )

save(list = c("ic.mean.P.matrices", "ic.mean.G.matrices", "all.lines.P.matrices", "all.lines.G.matrices"), file = "data/Matrices_dist_info.RData")

