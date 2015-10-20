source("./R/read_ratones.R")

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

lin_data = ldply(main.data, function(x) x$full) 
current.data <- select(lin_data, ID, .id, IS_PM:BA_OPI)

lm.within = lm(as.matrix(select(current.data, IS_PM:BA_OPI))~select(current.data, .id)[,1])


current.data_projected_W = cbind(select(current.data, .id, ID),  
                                 ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat))

Bmat = cov(daply(current.data_projected_W, .(.id), function(x) colMeans(x[,-c(1, 2)])))

resp <- cbind(select(current.data, .id, ID),  scale(ProjetaDados(select(current.data, IS_PM:BA_OPI), Wmat) %*% eigen(Bmat)$vectors[,1:3], scale = TRUE))
names(resp) <- c(".id", "ID", "CV1", "CV2", "CV3")
hulls <-ddply(resp, .(.id), plyr::summarise, "hpc1"=CV1[chull(CV1,CV2)],
                                             "hpc2"=CV2[chull(CV1,CV2)])
hulls %<>% separate(.id, c('treatment', 'line'))
resp %<>% separate(.id, c('treatment', 'line'))
cv_plot_12 <- ggplot(resp, aes(CV1, CV2)) +
  geom_polygon(aes(hpc1, hpc2, fill = line, group= interaction(line, treatment)), hulls, alpha=.3) + 
  geom_point(data = ddply(resp, .(line, treatment), numcolwise(mean)),
             aes(CV1, CV2, group= interaction(treatment, line), color = line, shape = treatment), size = 10) + 
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

save_plot(plot = cv_plot_12, "~/Dropbox/labbio/Shared Lab/Ratones_shared/figure1.pdf", base_height = 10)