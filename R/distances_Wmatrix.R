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
PRCOMP %>% biplot
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
  scale_fill_manual(values = c(c, h, s)) + scale_color_manual(values = c(c, h, s)) + theme_bw() + ggtitle("Cranial traits PC scores")
  
p49_plot <- lin_data %>% separate(.id, c('treatment', 'strain')) %>%
  ggplot(aes(treatment, P49, fill = strain)) + geom_boxplot() + scale_fill_manual(values = c(c, h, s)) + 
  facet_wrap(~SEX) + theme_bw() + ggtitle("Weigth at 49 days")
 
plot_grid(p49_plot, pc_plot, labels = c("A", "B"))
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/peso_pc.png")
plot_grid(global_stats_plot, pc_plot)
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/stats_pc.png")
plot_grid(DzPC1 , pc_plot)
ggsave("~/Dropbox/labbio/Shared Lab/Ratones_shared/dzpc1_pc.png")


