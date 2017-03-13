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

save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure2.pdf", figure_2, base_aspect_ratio = 1.3, base_height = 4.8)