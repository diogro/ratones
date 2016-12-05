if(!require(doMC)) {install.packages('doMC'); library(doMC)}
if(!require(viridis)) install.packages("viridis")
library(viridis)
registerDoMC(4)

#source("R/run_ratones_MCMCglmm_P.R")
source("R/read_ratones.R")

mcmc_stats = tbl_df(ldply(g_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE))
save(mcmc_stats, file = "./Rdatas/P_mcmc_stats")
load("./Rdatas/P_mcmc_stats")
names(mcmc_stats) <- gsub("pc1%", "pc1.percent", names(mcmc_stats))

  global_stats <- mcmc_stats %>% dplyr::select(.id, MeanSquaredCorrelation, flexibility, pc1.percent, evolvability) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
  global_stats$line <- gsub('control', 't', global_stats$line)
  global_stats$line <- factor(global_stats$line, levels = lines)
  r2_plot <- ggplot(filter(global_stats, variable == "MeanSquaredCorrelation"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(x = "Mean squared correlation", y = "Density") + theme(legend.position = c(0.8, 0.8), text = element_text(size = 20))
  flexibility_plot <- ggplot(filter(global_stats, variable == "flexibility"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line") + background_grid(major = 'y', minor = "none") +  panel_border() + labs(x = "Mean flexibility", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
  pc1.percent_plot <- ggplot(filter(global_stats, variable == "pc1.percent"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line") + background_grid(major = 'y', minor = "none") +  panel_border() + labs(x = "Proportion of variation\n in PC1", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
  evolvability_plot <- ggplot(filter(global_stats, variable == "evolvability"), aes(value, group = interaction(selection, variable, line), fill = interaction(selection, line))) + geom_density(alpha = 0.5) + scale_fill_manual(values = viridis(5), name = "Line")+ background_grid(major = 'y', minor = "none") +  panel_border() + labs(x = "Mean evolvability", y = "Density") + theme(legend.position = "none", text = element_text(size = 20))
figure_3 <- ggdraw() +
  draw_plot(r2_plot, 0, 0.5, 0.5, 0.5) +
  draw_plot(pc1.percent_plot, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(flexibility_plot, 0, 0, 0.5, 0.5) +
  draw_plot(evolvability_plot, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 20)

directionalVariation <- function(cov.matrix, line){
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(cov.matrix, ret.dim = 13)[[1]], delta_Z)
  data.frame(evolDZ = Evolvability(cov.matrix, Normalize(delta_Z)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             evolBeta = Evolvability(cov.matrix, Normalize(beta_NN)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])))
}

stats <- ldply(g_models, function(model) adply(model$Ps, 1,
                                               directionalVariation,
                                               model$line), .parallel = TRUE)

DzPC1 <- stats %>% dplyr::select(-X1) %>% ggplot(aes(DZpc1, group = .id, fill = .id)) + geom_density(alpha = 0.5) + labs(x = expression(paste("Vector correlation of ", Delta, "z and PC1")), y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'y', minor = "none") +  panel_border()+ theme(legend.position = c(0.15, 0.8), text = element_text(size = 20))

evolDZ <- stats %>% dplyr::select(-X1) %>% ggplot(aes(evolDZ, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  labs(x = "Scaled directional evolvability", y = "Density") + scale_fill_manual(values = viridis(5), name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) + background_grid(major = 'y', minor = "none") +  panel_border()+ theme(legend.position = "none", text = element_text(size = 20))

figure_4 <- ggdraw() +
  draw_plot(DzPC1 + theme(axis.title.x = element_text(size = rel(1.4))), 0, 0, 0.5, 1) +
  draw_plot(evolDZ+ theme(axis.title.x = element_text(size = rel(1.4))), 0.5, 0, 0.5, 1) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size = 20)

#evolBeta <- stats %>% dplyr::select(-X1) %>% ggplot(aes(evolBeta, group = .id, fill = .id)) + geom_density(alpha = 0.5) +  ggtitle(expression(paste("Ratio between mean evolvability and in the direction of ",Delta, "z"))) + labs(x = "Evolvability ratio") + background_grid(major = 'y', minor = "none") +  panel_border()

reps_RS = laply(g_models, function(x) mean(RandomSkewers(alply(x$Ps, 1), x$P)$correlation))
RS = llply(g_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps_RS)
rs_data <- RS[[1]]

reps_krz = laply(g_models, function(x) mean(KrzCor(alply(x$Ps, 1), x$P)$krz))
krz_data = llply(g_models, function(x) x$P) %>% KrzCor(repeat.vector = reps_krz)

library(xtable)
reps = rbind(reps_krz, reps_RS)
colnames(reps) <- c("Control t", "Downwards h", "Downwards s", "Upwards h'", "Upwards s'")
rownames(reps) <- c("Krzanowski", "Random Skewers")
xtable(t(reps), digits = 3)

mat_data <- rs_data
mat_data[lower.tri(mat_data)] <- t(krz_data)[lower.tri(krz_data)]
diag(mat_data) <- NA

myPalette <- viridis(50)
#myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
m.rs = melt(mat_data)
m.rs$Var1<- factor(m.rs$Var1, levels = levels(m.rs$Var1)[5:1])
m.rs.position = m.rs
m.rs.position$Var1 <- as.numeric(m.rs.position$Var1)
m.rs.position$Var2 <- as.numeric(m.rs.position$Var2)
m.rs.position$value= round(m.rs.position$value, 3)
m.rs.position$value[is.na(m.rs.position$value)] <- c("Control t", "Downwards h", "Downwards s", "Upwards h'", "Upwards s'")
figure_2 <- ggplot (m.rs) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette, limits = c(0.7, 1)) +
    ylab ('') + xlab ('') +
    geom_text(data = m.rs.position, size = 4, aes(x = Var2, y = Var1, label = value)) +
  theme_grey(base_size = 12, base_family = "") %+replace%
  theme(rect = element_rect(fill = "transparent", colour = NA,
                            color = NA, size = 0, linetype = 0),
        line = element_blank(),
        title = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = grid::unit(0, "lines"))

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure2.pdf", figure_2,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3, base_height = 5)

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure3.pdf", figure_3,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 5)

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure4.pdf", figure_4,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 5)

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure2.png", figure_2,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3, base_height = 5)

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure3.png", figure_3,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 5)

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure4.png", figure_4,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 5)

PCones <- t(laply(g_models, function(x) eigen(x$P)$vectors[,1]))
colnames(PCones) <- c("Control t", "Downwards h", "Downwards s", "Upwards h'", "Upwards s'")
rownames(PCones) <- colnames(main.data[[1]]$ed)

library(xtable)
xtable(PCones, digits = 3)

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

