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

#scaled_Ps = llply(names(r_models), function(group) aaply(r_models[[group]]$Ps, 1, function(x) x / (main.data[[group]]$ed.means %*% t(main.data[[group]]$ed.means))))
#names(scaled_Ps) <- names(r_models)
#scaled_mcmc_stats = tbl_df(ldply(scaled_Ps, function(x) adply(x, 1, MeanMatrixStatistics), .parallel = TRUE))
#x = cbind(mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability), scaled_mcmc_stats %>% select(evolvability))
#names(x)[5] <- "scaled_evolvability"
global_stats <- mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>% separate(.id, c( 'treatment', 'strain'))
#global_stats <- x %>% melt %>% separate(.id, c( 'treatment', 'strain'))
#{levels(global_stats$variable) <- c("Mean squared correlation", "Mean flexibility", "Mean evolvability", "Mean scaled evolvability")} 
levels(global_stats$variable) <- c("Mean squared correlation", "Mean flexibility", "Mean evolvability") 
  global_stats_plot <- ggplot(global_stats, aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() +  facet_wrap(~variable, scale = 'free') + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "", x = "Treatment") + ggtitle("Evolutionary statistics")

myPalette <- colorRampPalette(c("yellow", "white", "red"))(n = 100)
myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
m.rs = melt(mat_data) 
m.rs$Var1<- factor(m.rs$Var1, levels = levels(m.rs$Var1)[5:1])
m.rs.position = m.rs
m.rs.position$Var1 <- as.numeric(m.rs.position$Var1)
m.rs.position$Var2 <- as.numeric(m.rs.position$Var2)
m.rs.position$value= round(m.rs.position$value, 3)
m.rs.position$value[is.na(m.rs.position$value)] <- c("Control", "Increase h", "Increase s", "Reduce h", "Reduce s")
matrix_comparisons <- ggplot (m.rs) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette) +
    ylab ('') + xlab ('') + labs(title = "Matrix comparisons") + 
    geom_text(data = m.rs.position, size = 3, aes(x = Var2, y = Var1, label = value)) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_line(size = 0),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          rect = element_blank(), line = element_blank())
  
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
             evolDZ = Evolvability(cov.matrix, Normalize(delta_Z)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
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

DzPC1 <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'DZpc1') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() +  ggtitle(expression(paste("Correlation of ", Delta, "z and PC1"))) + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation") + background_grid(major = 'y', minor = "none") +  panel_border()

evolDZ <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'evolDZ') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() +  ggtitle(expression(paste("Ratio between mean evolvability and in the direction of ",Delta, "z"))) + scale_fill_manual(values = c(h, s)) + labs(y = "Evolvability ratio") + background_grid(major = 'y', minor = "none") +  panel_border()

figure_3 <- ggdraw() +
  draw_plot(global_stats_plot, 0, .5, 1, .5) +
  draw_plot(DzPC1, 0, 0, .5, .5) +
  draw_plot(matrix_comparisons, 0.5, 0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 20)
save_plot("~/Desktop/plot2by2.pdf", figure_3,
          ncol = 2, 
          nrow = 2, 
          base_aspect_ratio = 1.3
)

corDZDZ <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'corDZDZ') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0.32) +
  ggtitle("Correlation of observed change\n and expected change in multivariate mean") + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation")

treatment %<>%  separate(.id, c('treatment', 'strain'))
normDZ_DzPC1 = ggplot(treatment, aes(normDZ, DZpc1, group = interaction(treatment, strain), color = strain)) + geom_violin(aes(fill = strain), alpha = 0.3) + geom_jitter(aes(shape = treatment), size = 3, position = position_jitter(width = .03)) + scale_fill_manual(values = c(h, s)) + scale_color_manual(values = c(h, s))




