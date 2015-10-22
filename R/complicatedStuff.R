traits = select(full_data, P49, IS_PM:BA_OPI)
traits$P49 %<>% log
#delta_Zs <- llply(main.data, function(x) x$ed.means - main.data$control.control$ed.means)
delta_Zs <- llply(main.data, function(x) x$ed.means - colMeans(traits))
plsr <- llply(main.data, function(x) x$plsr)

directionalVariation <- function(cov.matrix, line){
  delta_Z <- delta_Zs[[line]]
  beta_s <- solve(cov.matrix, delta_Z)
  beta_NN <- solve(ExtendMatrix(cov.matrix, ret.dim = 13)[[1]], delta_Z)  
  beta = c(1, rep(0, 35))
  data.frame(corDZDZ = abs(vectorCor(cov.matrix %*% beta, delta_Z)),
             evolDZ = Evolvability(cov.matrix, Normalize(delta_Z)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])),
             normDZ = Norm(delta_Z),
             normPLS = Norm(main.data[[line]]$plsr))
}

treatment <- ldply(r_models[-1], function(model) adply(model$Ps, 1, 
                                                       directionalVariation, 
                                                       model$line), .parallel = TRUE)
treatment$type <- 'treatment'
control   <- ldply(r_models[-1], function(model) adply(r_models[['control.control']]$Ps, 1, 
                                                       directionalVariation, 
                                                       model$line), .parallel = TRUE)
control$type <- 'control'
stats <- melt(rbind(treatment, control))[-2]

DzPC1 <- stats %>% separate(.id, c( 'treatment', 'line')) %>% filter(variable == 'DZpc1') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, line, type), fill = line)) + geom_boxplot() + scale_fill_manual(values = c(h, s)) + labs(y = expression(paste("Vector correlation of ", Delta, "z and PC1"))) + background_grid(major = 'y', minor = "none") +  panel_border()

evolDZ <- stats %>% separate(.id, c( 'treatment', 'line')) %>% filter(variable == 'evolDZ') %>% filter(type == "treatment") %>%   ggplot(aes(treatment, value, group = interaction(treatment, line, type), fill = line)) + geom_boxplot() +  ggtitle(expression(paste("Ratio between mean evolvability and in the direction of ",Delta, "z"))) + scale_fill_manual(values = c(h, s)) + labs(y = "Evolvability ratio") + background_grid(major = 'y', minor = "none") +  panel_border()


corDZDZ <- stats %>% separate(.id, c( 'treatment', 'line')) %>% filter(variable == 'corDZDZ') %>% filter(type == "treatment") %>% 
  ggplot(aes(treatment, value, group = interaction(treatment, line, type), fill = line)) + geom_boxplot() + theme_bw() + geom_hline(yintercept = 0.32) +
  ggtitle("Correlation of observed change\n and expected change in multivariate mean") + scale_fill_manual(values = c(h, s)) + labs(y = "Vector correlation")

treatment %<>%  separate(.id, c('treatment', 'line'))

normDZ_DzPC1 = ggplot(treatment, aes(normDZ, DZpc1, group = interaction(treatment, line), color = line)) + geom_violin(aes(fill = line), alpha = 0.3) + geom_jitter(aes(shape = treatment), size = 3, position = position_jitter(width = .03)) + scale_fill_manual(values = c(h, s)) + scale_color_manual(values = c(h, s))

stats %>% separate(.id, c( 'treatment', 'line')) %>% filter(variable == 'normDZ' | variable == 'normPLS') %>% filter(type == "treatment") %>% ddply(.(variable, line, treatment), numcolwise(mean)) %>% spread(variable, value) %>%
  ggplot(aes(normDZ, normPLS, group = interaction(treatment, line), color = interaction(treatment, line))) + geom_point() 


normalizedEvolvability <- function(cov.matrix, line){
  delta_Z <- delta_Zs[[line]]
  cov.matrix = cov.matrix / main.data[[line]]$gm_mean
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
                                                       model$line), .parallel = TRUE)
treatment$type <- 'treatment'
control   <- ldply(r_models[-1], function(model) adply(r_models[['control.control']]$Ps, 1, 
                                                       normalizedEvolvability, 
                                                       model$line), .parallel = TRUE)
control$type <- 'control'
evol_norm <- melt(rbind(treatment, control))[-2]


evol_norm %>% separate(.id, c( 'treatment', 'line')) %>% filter(variable == 'evol_Random_notDZ' | variable == 'evol_DZ') %>%
  ggplot(aes(type, value, group = interaction(treatment, line, variable, type), fill = line)) + geom_boxplot() +  facet_grid(variable~treatment, scales = 'free_y') + scale_fill_manual(values = c(h, s)) + background_grid(major = 'y', minor = "none") +  panel_border() + labs(y = "", x = "Treatment")