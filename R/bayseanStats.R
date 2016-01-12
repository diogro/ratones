if(!require(doMC)) {install.packages('doMC'); library(doMC)}
registerDoMC(4)

matrix_median <- function(mat.list){
  x = laply(mat.list, identity)
  aaply(x, 2:3, median)
}

#source("R/run_ratones_MCMCglmm_P.R")
source("R/read_ratones.R")

r_models = llply(main.data, function(x) evolqg:::CalculateMatrix_Baysean(x$model, samples = 100, nu = 2))
for(i in 1:length(r_models)){
  r_models[[i]]$line <- names(r_models)[1]
  r_models[[i]]$P <- matrix_median(r_models[[i]]$Ps)
}

mcmc_stats = tbl_df(ldply(r_models, function(x) ldply(x$Ps, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE, .inform = TRUE))
# mcmc_stats = tbl_df(ldply(r_models, function(x) adply(x$Ps, 1, MeanMatrixStatistics, .progress = "text"), .parallel = TRUE, .inform = TRUE))
# save(mcmc_stats, file = "./Rdatas/mcmc_stats")
#load("./Rdatas/mcmc_stats")
names(mcmc_stats) <- gsub("pc1%", "pc1.percent", names(mcmc_stats))

  global_stats <- mcmc_stats %>% select(.id, MeanSquaredCorrelation) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
  global_stats$line <- gsub('control', 't', global_stats$line)
  global_stats$line <- factor(global_stats$line, levels = lines)
  r2_plot <- ggplot(global_stats, aes(line, value, group = interaction(selection, line, variable), fill = selection)) + geom_boxplot() + scale_fill_manual(values = c(c, dw, up)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean squared correlation", x = "") + theme(legend.position = c(0.15, 0.8))

  global_stats <- mcmc_stats %>% select(.id, flexibility) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
  global_stats$line <- gsub('control', 't', global_stats$line)
  global_stats$line <- factor(global_stats$line, levels = lines)
  flexibility_plot <- ggplot(global_stats, aes(line, value, group = interaction(selection, line, variable), fill = selection)) + geom_boxplot() + scale_fill_manual(values = c(c, dw, up)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean flexibility", x = "") + theme(legend.position = "none", text = element_text(size = 20))
  
  global_stats <- mcmc_stats %>% select(.id, pc1.percent) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
  global_stats$line <- gsub('control', 't', global_stats$line)
  global_stats$line <- factor(global_stats$line, levels = lines)
  pc1.percent_plot <- ggplot(global_stats, aes(line, value, group = interaction(selection, line, variable), fill = selection)) + geom_boxplot() + scale_fill_manual(values = c(c, dw, up)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Proportion of variation\n in PC1", x = "") + theme(legend.position = "none", text = element_text(size = 20))


# reps_RS = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = RandomSkewers)
# RS = llply(r_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps_RS)
# rs_data <- RS[[1]]
# 
# reps_krz = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = KrzCor)
# krz_data = llply(r_models, function(x) x$P) %>% KrzCor(repeat.vector = reps_krz)
# 
# library(xtable)
# reps = rbind(reps_krz, reps_RS)
# colnames(reps) <- c("Control", "Increase h", "Increase s", "Reduce h", "Reduce s")
# rownames(reps) <- c("Krzanowski", "Random Skewers")
# xtable(reps, digits = 3)
# 
# mat_data <- rs_data
# mat_data[lower.tri(mat_data)] <- t(krz_data)[lower.tri(krz_data)]
# diag(mat_data) <- NA
# save(mat_data,file = "./Rdatas/mat_data.Rdata")
#load("./Rdatas/mat_data.Rdata")

reps_RS = laply(r_models, function(x) mean(RandomSkewers(x$Ps, x$MAP)[,1]))
RS = llply(r_models, function(x) x$MAP) %>% RandomSkewers(repeat.vector = reps_RS)
rs_data <- RS[[1]]
  
reps_krz = laply(r_models, function(x) mean(KrzCor(x$Ps, x$MAP)))
krz_data = llply(r_models, function(x) x$MAP) %>% KrzCor(repeat.vector = reps_krz)
  
library(xtable)
reps = rbind(reps_krz, reps_RS)
colnames(reps) <- c("Control t", "Increase h", "Increase s", "Reduce h", "Reduce s")
rownames(reps) <- c("Krzanowski", "Random Skewers")
xtable(reps, digits = 3)
  
mat_data <- rs_data
mat_data[lower.tri(mat_data)] <- t(krz_data)[lower.tri(krz_data)]
diag(mat_data) <- NA
  
myPalette <- colorRampPalette(c("yellow", "white", "red"))(n = 100)
#myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
m.rs = melt(mat_data)
m.rs$Var1<- factor(m.rs$Var1, levels = levels(m.rs$Var1)[5:1])
m.rs.position = m.rs
m.rs.position$Var1 <- as.numeric(m.rs.position$Var1)
m.rs.position$Var2 <- as.numeric(m.rs.position$Var2)
m.rs.position$value= round(m.rs.position$value, 3)
m.rs.position$value[is.na(m.rs.position$value)] <- c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")
matrix_comparisons <- ggplot (m.rs) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette) +
    ylab ('') + xlab ('') +
    geom_text(data = m.rs.position, size = 4, aes(x = Var2, y = Var1, label = value)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_line(size = 0),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          rect = element_blank(), line = element_blank())

traits = select(full_data, IS_PM:BA_OPI)
#delta_Zs <- llply(main.data, function(x) x$ed.means - main.data$control.control$ed.means)
delta_Zs <- llply(main.data, function(x) x$ed.means - colMeans(traits))

directionalVariation <- function(cov.matrix, line){
  delta_Z <- delta_Zs[[line]]
  data.frame(DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])))
}


DzPC1_stat <- ldply(r_models[-1], function(model) ldply(model$Ps,
                                                          directionalVariation,
                                                          model$line), .parallel = TRUE)

DzPC1_data <- DzPC1_stat %>% select(.id, DZpc1) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
DzPC1_data$line = factor(DzPC1_data$line, levels = lines[-1])
DzPC1 = ggplot(DzPC1_data, aes(line, value, group = interaction(selection, line), fill = selection)) + geom_boxplot() + scale_fill_manual(values = c(dw, up)) + labs(y = expression(paste("Vector correlation of ", Delta, "z and PC1"))) + background_grid(major = 'y', minor = "none") +  panel_border() + theme_bw()

scaledEvolvability <- function(cov.matrix, line){
  delta_Z <- delta_Zs[[line]]
  cov.matrix = cov.matrix / main.data[[line]]$gm_mean

  data.frame(scaled_evol = mean(Evolvability(cov.matrix)))
}

scaled_mean_evol <- ldply(r_models, function(model) adply(model$Ps, 1,
                                                       scaledEvolvability,
                                                       model$line), .parallel = TRUE)

evolvability_data <- scaled_mean_evol %>% select(.id, scaled_evol) %>% melt %>% separate(.id, c('selection', 'line'), sep = "\\.")
evolvability_data$line <- gsub('control', 't', evolvability_data$line)
evolvability_data$line = factor(evolvability_data$line, levels = lines)
evolvability_plot <- ggplot(evolvability_data, aes(line, value, group = interaction(selection, line, variable), fill = selection)) + geom_boxplot() + scale_fill_manual(values = c(c, dw, up)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean evolvability", x = "") + theme(legend.position = "none", text = element_text(size = 20))


figure_2 <- ggdraw() +
  draw_plot(r2_plot + theme(axis.title.y = element_text(size = rel(1.4))), 0, 0.5, 0.5, 0.5) +
  draw_plot(pc1.percent_plot, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(flexibility_plot, 0, 0, 0.5, 0.5) +
  draw_plot(evolvability_plot, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 20)
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure2.pdf", figure_2,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3, base_height = 5)

figure_3 <- matrix_comparisons
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure3.pdf", figure_3,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 5)

figure_4 <- DzPC1 + theme(legend.position = c(0.2, 0.15))
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure4.pdf", figure_4,
          ncol = 1,
          nrow = 1,
          base_aspect_ratio = 1.3, base_height = 4)


PCones <- t(laply(r_models, function(x) eigen(x$MAP)$vectors[,1]))
colnames(PCones) <- c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")
rownames(PCones) <- rownames(main.data[[1]]$cov.matrix)

library(xtable)
xtable(PCones, digits = 3)

PC1_iso_cor = aaply(PCones, 2, function(x) abs(vectorCor(x, sqrt(rep(1/35, 35)))))

load("./Rdatas/PG.Calomys.RData")

G_comp = cbind(RandomSkewers(llply(r_models, `[[`, "MAP"), Calomys.Pacote$G)[,1:2],
                      KrzCor(llply(r_models, `[[`, "MAP"), Calomys.Pacote$G)[,2])
names(G_comp) = c(".id", "Random Skewers", "Krzanowski")
G_comp %>% separate(.id, c("selection", "line"), sep = "\\.") %>% {xtable(.[,c(2, 1, 3, 4)])}
 