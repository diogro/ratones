source("./R/bayseanStats.R")

#m.data = gather(raw.data, trait, value, IS_PM:BA_OPI)
#lmer(value ~ trait:SEX:LIN + (0 + trait|ID), m.data)
lin_data = ldply(main.data, function(x) x$full) 
individual_ID <- unlist(llply(main.data, function(x) x$info$ID))
myformula = paste0('~', paste(names(select(lin_data, IS_PM:BA_OPI)), collapse = '+'))
rownames(lin_data) <- individual_ID
#lin_data = filter(lin_data, .id != 'reduce.s', .id != 'increase.s') 
PRCOMP <- princomp(as.formula(myformula), data = lin_data, scale = T)
PRCOMP$scores
######## Biplot PC1 x PC2 ##############
current.data <- select(lin_data, ID, .id, IS_PM:BA_OPI)
rownames(current.data) <- rownames(lin_data)
#PRCOMP %>% biplot
#PRCOMP <- princomp(data.frame(na.omit(current.data$sizeless ) ) )
#resp <- current.data[rownames(PRCOMP$scores), ] #respectivos dados aos participantes da PCA
#resp %<>% mutate(., PC1 = PRCOMP$scores[,1], PC2=PRCOMP$scores[,2]) 
resp <- cbind(select(current.data, .id, ID),  as.matrix(select(current.data, IS_PM:BA_OPI)) %*% eigen(Wmat)$vectors[,1:2])
names(resp) <- c(".id", "ID", "PC1", "PC2")
hulls <-ddply(resp, .(.id), plyr::summarise, "hpc1"=PC1[chull(PC1,PC2)],
                                             "hpc2"=PC2[chull(PC1,PC2)])
hulls %<>% separate(.id, c('treatment', 'strain'))
resp %<>% separate(.id, c('treatment', 'strain'))
pc_plot <- ggplot(resp, aes(PC1, PC2)) +
  #geom_point(aes(PC1, PC2, shape = treatment), size = 4, alpha = 0.5) +
  geom_polygon(aes(hpc1, hpc2, fill = strain, group= interaction(strain, treatment)), hulls, alpha=.3) + 
  geom_point(data = ddply(resp, .(strain, treatment), numcolwise(mean)),
             aes(PC1, PC2, group= interaction(treatment, strain), color = strain, shape = treatment), size = 10) + 
  scale_fill_manual(values = c(c, h, s)) + scale_color_manual(values = c(c, h, s)) + theme_bw() + ggtitle("Cranial traits Within-group PC scores")
  
lin_data$SEX <- factor(lin_data$SEX)
levels(lin_data$SEX) <- c("Females", "Males")
p49_plot <- lin_data %>% separate(.id, c('treatment', 'strain')) %>%
  ggplot(aes(treatment, P49, fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + 
  facet_wrap(~SEX) + theme_bw() + ggtitle("Weigth at 49 days") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 30), 
        plot.title = element_text(size = 40)) 

global_stats_plot = global_stats_plot +   theme(text = element_text(size = 30),
                                                legend.text = element_text(size = 30), 
                                                plot.title = element_text(size = 40)) 


## Canonical variates
ProjetaDados = function(y,var.y){
  y = as.matrix(y)
  n = nrow(y)
  p = ncol(y)
  eigen.y = eigen(var.y)
  eVal = eigen.y$values
  eVec = eigen.y$vectors
  Scores = array(0., c(n,p))
  for( i in 1:n){
    Scores[i,] = t(eVec)%*%(as.numeric(y[i,]))
  }
  for(i in 1:p){
    Scores[,i] <- Scores[,i]/eVal[i]
  }
  return(Scores)
}

lm.within = lm(as.matrix(select(current.data, IS_PM:BA_OPI))~select(current.data, .id)[,1])

current.data_projected_W = cbind(select(current.data, .id, ID),  
                                 ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat))

Bmat = cov(daply(current.data_projected_W, .(.id), function(x) colMeans(x[,-c(1, 2)])))

resp <- cbind(select(current.data, .id, ID),  scale(ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat) %*% eigen(Bmat)$vectors[,1:3], scale = TRUE))
names(resp) <- c(".id", "ID", "CV1", "CV2", "CV3")
hulls <-ddply(resp, .(.id), plyr::summarise, "hpc1"=CV1[chull(CV1,CV2)],
                                             "hpc2"=CV2[chull(CV1,CV2)])
hulls %<>% separate(.id, c('treatment', 'strain'))
resp %<>% separate(.id, c('treatment', 'strain'))
cv_plot_12 <- ggplot(resp, aes(CV1, CV2)) +
  geom_polygon(aes(hpc1, hpc2, fill = strain, group= interaction(strain, treatment)), hulls, alpha=.3) + 
  geom_point(data = ddply(resp, .(strain, treatment), numcolwise(mean)),
             aes(CV1, CV2, group= interaction(treatment, strain), color = strain, shape = treatment), size = 10) + 
  scale_fill_manual(values = c(c, h, s)) + scale_color_manual(values = c(c, h, s)) + 
  theme(legend.position = "none", axis.text = element_text(size = 30), 
        axis.title = element_text(size = 30)) + 
  labs(x = "First canonical variate", y = "Second canonial variate") +
  annotate("text", 1.5, 1.5, label = "Increase", color = h, size = 15) + 
  annotate("text", -1.5, 1.5, label = "Increase", color = s, size = 15) + 
  annotate("text", 1, -2, label = "Reduce", color = h, size = 15) + 
  annotate("text", -1.5, -2, label = "Reduce", color = s, size = 15) + 
  annotate("text", 0.5, 2.5, label = "Control", color = c, size = 15) + 
  annotate("text", -2.25, 0.5, label = "s", color = s, size = 15) + 
  annotate("text", 1.5, 0, label = "h", color = h, size = 15) 
# 
# current.data_projected_W = cbind(select(current.data, .id, ID),  ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat))
# Bmat = cov(daply(current.data_projected_W, .(.id), function(x) colMeans(x[,-c(1, 2)])))
# resp <- cbind(select(current.data, .id, ID),  ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat) %*% eigen(Bmat)$vectors[,1:3])
# names(resp) <- c(".id", "ID", "CV1", "CV2", "CV3")
# hulls <-ddply(resp, .(.id), plyr::summarise, "hpc2"=CV2[chull(CV2,CV3)],
#                                              "hpc3"=CV3[chull(CV2,CV3)])
# resp %<>% separate(.id, c('treatment', 'strain'))
# hulls %<>% separate(.id, c('treatment', 'strain'))
# cv_plot_23<- ggplot(resp, aes(CV2, CV3)) +
#   #geom_point(aes(PC1, PC3, shape = treatment), size = 4, alpha = 0.5) +
#   geom_polygon(aes(hpc2, hpc3, fill = strain, group= interaction(strain, treatment)), hulls, alpha=.3) +
#   geom_point(data = ddply(resp, .(strain, treatment), numcolwise(mean)),
#              aes(CV2, CV3, group= interaction(treatment, strain), color = strain, shape = treatment), size = 10) + 
#   scale_fill_manual(values = c(c, h, s)) + scale_color_manual(values = c(c, h, s)) + theme_bw() + ggtitle("Cranial traits Canonical Variates Scores")


plot_grid(p49_plot, cv_plot_12, labels = c("A", "B"))
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/peso_pc.png", width = 15, height = 10)
plot_grid(global_stats_plot, cv_plot_12)
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/stats_pc.png", width = 15, height = 10)
plot_grid(DzPC1 , cv_plot_12)
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/dzpc1_pc.png", width = 15, height = 10)

ggsave(plot = p49_plot, "./md/p49_plot.png", width = 20, height = 15)
save_plot(plot = cv_plot_12, "./md/cv_plot_12.pdf", base_height = 10)
ggsave(plot = global_stats_plot, "./md/stats.png", width = 20, height = 15)

reps_RS = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = RandomSkewers)
RS = llply(r_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps_RS)
rs_data <- RS[[1]]

reps_krz = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = KrzCor)
krz_data = llply(r_models, function(x) x$P) %>% KrzCor(repeat.vector = reps_krz)

mat_data <- rs_data
mat_data[lower.tri(mat_data)] <- t(krz_data)[lower.tri(krz_data)]
diag(mat_data) <- NA
