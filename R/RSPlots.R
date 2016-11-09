library(evolqg)
library(magrittr)
library(ggplot2)
library(reshape2)
library(cowplot)

num_samples = 500

mats = llply(r_models, function(x) x$Ps)
maps = llply(r_models, function(x) x$MAP)

p <- length(mats)
n <- names(mats)

RSresults <- list()
for(i in 1:(p-1)){
  for(j in (i+1):p){
    pair <- paste(n[i], n[j], sep = "-")
    RSresults[[pair]] <- RandomSkewers(mats[[i]], mats[[j]], parallel = TRUE)
  }
}
rs_intervals <- ldply(RSresults, quantile, c(0.025, 0.5, 0.975))
names(rs_intervals) <- c("pair", "lower", "median", "upper")
rs_intervals$comparison <- "Random Skewers"

rs.map = melt(RandomSkewers(maps)[1]) %>% filter(value > 0)
rs.map$pair = paste(rs.map[,2], rs.map[,1], sep = "-")
rs.map$comparison = "Random Skewers"
rs.map %<>% select(value, pair, comparison)

KRZresults <- list()
for(i in 1:(p-1)){
  for(j in (i+1):p){
    pair <- paste(n[i], n[j], sep = "-")
    KRZresults[[pair]] <- KrzCor(mats[[i]], mats[[j]])
  }
}
krz_intervals <- ldply(KRZresults, quantile, c(0.025, 0.5, 0.975))
names(krz_intervals) <- c("pair", "lower", "median", "upper")
krz_intervals$comparison <- "Krzanowski"

krz.map = melt(KrzCor(maps)) %>% filter(value > 0)
krz.map$pair = paste(krz.map[,2], krz.map[,1], sep = "-")
krz.map$comparison = "Krzanowski"
krz.map %<>% select(value, pair, comparison)

intervals <- rbind(rs_intervals, krz_intervals)
comp.map <- rbind(rs.map, krz.map)

mean_comp <- ddply(intervals, .(comparison), numcolwise(mean))

rs_plot <- ggplot(intervals, aes_string("pair", "median", group = "comparison", color = "comparison")) + 
  geom_point(position = position_dodge(width = 0.4), aes_string(shape = "comparison"), size = 3) + 
  geom_point(data = comp.map, position = position_dodge(width = 0.4), aes_string(y = "value", shape = "comparison"), size = 3) + 
  geom_linerange(position = position_dodge(width = 0.4), 
                 aes_string(ymin = "lower", ymax = "upper", position = "comparison")) +  
  geom_hline(data = mean_comp, aes_string(yintercept = "median", color = "comparison")) + 
  coord_flip() + theme(legend.position = c(0.2, 0.9)) + labs(y = "Matrix comparison", x = "Matrix pair")
save_plot("~/Desktop/comparison.png", rs_plot, base_height = 5, base_aspect_ratio = 2)
