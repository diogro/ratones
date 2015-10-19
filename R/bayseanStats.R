source("R/run_ratones_MCMCglmm_P.R")


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

# global_stats <- mcmc_stats %>% select(.id, MeanSquaredCorrelation, flexibility, evolvability) %>% melt %>% separate(.id, c( 'treatment', 'strain'))
# levels(global_stats$variable) <- c("Mean squared correlation", "Mean flexibility", "Mean evolvability") 
#   global_stats_plot <- ggplot(global_stats, aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() +  facet_wrap(~variable, scale = 'free') + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "", x = "") + ggtitle("Evolutionary statistics")

  global_stats <- mcmc_stats %>% select(.id, MeanSquaredCorrelation) %>% melt %>% separate(.id, c('treatment', 'strain'))
  r2_plot <- ggplot(global_stats, aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean squared correlation", x = "") 
  
  global_stats <- mcmc_stats %>% select(.id, flexibility) %>% melt %>% separate(.id, c('treatment', 'strain'))
  flexibility_plot <- ggplot(global_stats, aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean flexibility", x = "") + theme(legend.position = "none", text = element_text(size = 20))
  
  global_stats <- mcmc_stats %>% select(.id, pc1.percent) %>% melt %>% separate(.id, c('treatment', 'strain'))
  pc1.percent_plot <- ggplot(global_stats, aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Proportion of variation\n in PC1", x = "") + theme(legend.position = "none", text = element_text(size = 20))
  
  
# reps_RS = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = RandomSkewers)
# RS = llply(r_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps_RS)
# rs_data <- RS[[1]]
# 
# reps_krz = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = KrzCor)
# krz_data = llply(r_models, function(x) x$P) %>% KrzCor(repeat.vector = reps_krz)
# 
# mat_data <- rs_data
# mat_data[lower.tri(mat_data)] <- t(krz_data)[lower.tri(krz_data)]
# diag(mat_data) <- NA
# save(mat_data,file = "./Rdatas/mat_data.Rdata")
load("./Rdatas/mat_data.Rdata")

myPalette <- colorRampPalette(c("yellow", "white", "red"))(n = 100)
#myPalette <- colorRampPalette(c("yellow", "white", "purple"))(n = 100)
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
    ylab ('') + xlab ('') +
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
#delta_Zs <- llply(main.data, function(x) x$ed.means - main.data$control.control$ed.means)
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
             normDZ = Norm(delta_Z),
             normPLS = Norm(main.data[[strain]]$plsr))
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

DzPC1 <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'DZpc1') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(h, s)) + labs(y = expression(paste("Vector correlation of ", Delta, "z and PC1"))) + background_grid(major = 'y', minor = "none") +  panel_border()

evolDZ <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'evolDZ') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() +  ggtitle(expression(paste("Ratio between mean evolvability and in the direction of ",Delta, "z"))) + scale_fill_manual(values = c(h, s)) + labs(y = "Evolvability ratio") + background_grid(major = 'y', minor = "none") +  panel_border()


corDZDZ <- stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'corDZDZ') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, strain, type), fill = strain)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0.32) +
  ggtitle("Correlation of observed change\n and expected change in multivariate mean") + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation")

treatment %<>%  separate(.id, c('treatment', 'strain'))

normDZ_DzPC1 = ggplot(treatment, aes(normDZ, DZpc1, group = interaction(treatment, strain), color = strain)) + geom_violin(aes(fill = strain), alpha = 0.3) + geom_jitter(aes(shape = treatment), size = 3, position = position_jitter(width = .03)) + scale_fill_manual(values = c(h, s)) + scale_color_manual(values = c(h, s))

stats %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'normDZ' | variable == 'normPLS') %>% filter(type == "treatment") %>% ddply(.(variable, strain, treatment), numcolwise(mean)) %>% spread(variable, value) %>%
  ggplot(aes(normDZ, normPLS, group = interaction(treatment, strain), color = interaction(treatment, strain))) + geom_point() 


normalizedEvolvability <- function(cov.matrix, strain){
  delta_Z <- delta_Zs[[strain]]
  cov.matrix = cov.matrix / main.data[[strain]]$gm_mean
  pc1 = eigen(cov.matrix)$vectors[,1]
  beta_mat_notDZ = array(NA, dim = c(length(delta_Z), num_vector <- 1000))
  for(i in 1:num_vector){
    while(TRUE){
      rand_vector = rnorm(length(delta_Z))
      if(abs(vectorCor(rand_vector, delta_Z)) < 0.01){
        beta_mat_notDZ[, i] <- Normalize(rand_vector)
        break
      }
    }
  }
  beta_mat_notPC1 = array(NA, dim = c(length(delta_Z), num_vector <- 1000))
  for(i in 1:num_vector){
    while(TRUE){
      rand_vector = rnorm(length(delta_Z))
      if(abs(vectorCor(rand_vector, pc1)) < 0.01){
        beta_mat_notPC1[, i] <- Normalize(rand_vector)
        break
      }
    }
  }
  data.frame(evol_Random_notDZ = mean(Evolvability(cov.matrix, beta_mat_notDZ)),
             evol_DZ = Evolvability(cov.matrix, Normalize(delta_Z)),
             evol_Random_notPC1 = mean(Evolvability(cov.matrix, beta_mat_notPC1)),
             evol_PC1 = Evolvability(cov.matrix, pc1), 
             scaled_evol = mean(Evolvability(cov.matrix)))
}
treatment <- ldply(r_models[-1], function(model) adply(model$Ps, 1, 
                                                       normalizedEvolvability, 
                                                       model$strain), .parallel = TRUE)
treatment$type <- 'treatment'
control   <- ldply(r_models[-1], function(model) adply(r_models[['control.control']]$Ps, 1, 
                                                       normalizedEvolvability, 
                                                       model$strain), .parallel = TRUE)
control$type <- 'control'
evol_norm <- melt(rbind(treatment, control))[-2]


evol_norm %>% separate(.id, c( 'treatment', 'strain')) %>% filter(variable == 'evol_Random_notDZ' | variable == 'evol_DZ') %>%
ggplot(aes(type, value, group = interaction(treatment, strain, variable, type), fill = strain)) + geom_boxplot() +  facet_grid(variable~treatment, scales = 'free_y') + scale_fill_manual(values = c(h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "", x = "Treatment")

scaledEvolvability <- function(cov.matrix, strain){
  delta_Z <- delta_Zs[[strain]]
  cov.matrix = cov.matrix / main.data[[strain]]$gm_mean

  data.frame(scaled_evol = mean(Evolvability(cov.matrix)))
}

scaled_mean_evol <- ldply(r_models, function(model) adply(model$Ps, 1, 
                                                       scaledEvolvability, 
                                                       model$strain), .parallel = TRUE)

evolvability_plot <- scaled_mean_evol %>% select(.id, scaled_evol) %>% melt %>% separate(.id, c('treatment', 'strain')) %>% ggplot(aes(treatment, value, group = interaction(treatment, strain, variable), fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "Mean evolvability", x = "") + theme(legend.position = "none", text = element_text(size = 20))


figure_2 <- ggdraw() +
  draw_plot(r2_plot + 
              theme(legend.position = c(0.15, 0.8)) , 0, 0.5, 0.5, 0.5) +
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

figure_4 <- DzPC1 + theme(legend.position = c(0.1, 0.91))
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure4.pdf", figure_4,
          ncol = 1, 
          nrow = 1, 
          base_aspect_ratio = 1.3, base_height = 5)

