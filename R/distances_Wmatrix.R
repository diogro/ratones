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
  scale_fill_manual(values = c(c, h, s)) + scale_color_manual(values = c(c, h, s)) + theme_bw() + 
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 30), 
        plot.title = element_text(size = 30)) + 
  annotate("text", 1.5, 1.5, label = "Increase", color = h, size = 20) + 
  annotate("text", -1.5, 1.5, label = "Increase", color = s, size = 20) + 
  annotate("text", 1, -2, label = "Reduce", color = h, size = 20) + 
  annotate("text", -1.5, -2, label = "Reduce", color = s, size = 20) + 
  annotate("text", 0.5, 2.5, label = "Control", color = c, size = 20) + 
  annotate("text", -2.25, 0.5, label = "s", color = s, size = 20) + 
  annotate("text", 1.5, 0, label = "h", color = h, size = 20) + 
  ggtitle("Cranial traits Canonical Variates Scores")
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
ggsave(plot = cv_plot_12, "./md/cv_plot_12.png", width = 20, height = 15)
ggsave(plot = global_stats_plot, "./md/stats.png", width = 20, height = 15)

reps = llply(r_models, function(x) x$P) %>% laply(., MonteCarloRep, sample.size = 50, ComparisonFunc = RandomSkewers)
RS = llply(r_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps)
mat_data <- RS[[1]]

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}



#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0.6, 0.8, length=33),  # for red
               seq(0.8, 0.9, length=34),
               seq(0.9,   1,  length=34))              # for green

#creates a 5 x 5 inch image
png("./md/heatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(mat_data, 
          cellnote = round(mat_data, 3),  # same data set for cell labels
          #main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          #dendrogram="row",     # only draw a row dendrogram
          Colv="NA", Rowv = "NA")            # turn off column clustering

dev.off()               # close the PNG device


llply(r_models, function(x) x$P) %>% RandomSkewers(repeat.vector = reps) -> RS
color2D.matplot(RS[[1]])

