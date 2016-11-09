MeanP49<- read.csv("mean_P49.csv", header=TRUE, as.is = T)
MeanP49$GEN <- as.factor(MeanP49$GEN)

MeanP49 %<>% melt(value.name = "P49") 

MeanP49$selection <- MeanP49$variable

MeanP49$selection %<>% gsub("s", "downwards", ., ignore.case = F)
MeanP49$selection %<>% gsub("t", "control", ., ignore.case = F)
MeanP49$selection %<>% gsub("h", "downwards", ., ignore.case = F)
MeanP49$selection %<>% gsub("downwards.", "upwards", .)
unique(MeanP49$selection)

names(MeanP49)[3] <- "line"
MeanP49$original.line <- MeanP49$line
MeanP49$original.line %<>% gsub("\\.", "'", .)


hand.plot <- MeanP49 %>% filter(GEN == 55 & original.line == "t" |
                   GEN == 47 & original.line == "h'" |
                   GEN == 47 & original.line == "h" |
                   GEN == 48 & original.line == "h" |
                   GEN == 52 & original.line == "s'" |
                   GEN == 52 & original.line == "s" |
                   GEN == 53 & original.line == "s" 
                   )

Demog<- read.csv("Demografia.csv", header=TRUE)
Demog$GEN <- as.factor(Demog$GEN)
Demog %<>% melt(value.name = "Census")
Demog$selection <- Demog$variable
Demog$selection %<>% gsub("s", "downwards", ., ignore.case = F)
Demog$selection %<>% gsub("t", "control", ., ignore.case = F)
Demog$selection %<>% gsub("h", "downwards", ., ignore.case = F)
Demog$selection %<>% gsub("downwards.", "upwards", ., ignore.case = F)
unique(Demog$selection)

names(Demog)[3] <- "line"
Demog$original.line <- Demog$line
Demog$original.line %<>% gsub("\\.", "'", .)

hand.plot.dem <- Demog %>% filter(GEN == 55 & original.line == "t" |
                                  GEN == 47 & original.line == "h'" |
                                  GEN == 47 & original.line == "h" |
                                  GEN == 48 & original.line == "h" |
                                  GEN == 52 & original.line == "s'" |
                                  GEN == 52 & original.line == "s" |
                                  GEN == 53 & original.line == "s" 
)

labels.lines = c("t" = "control t", "h" = "downwards h", "h'" = "upwards h'", "s" = "downwards s", "s'" = "upwards s'")

demo.plot <- Demog %>% ggplot(aes(x= GEN, y= Census, group = interaction(SEX, line), color = SEX, shape =original.line) ) +
  geom_point(size =3) +
  geom_line() +
  geom_point(data = hand.plot.dem, shape = 16, color = "darkgrey", size = 5, alpha = 0.8) +
  facet_wrap(~ original.line, nrow = 3, scales =  "free_y", labeller = as_labeller(labels.lines) )+
  scale_x_discrete(breaks = seq(from = 0, to = 55, 5), labels = seq(from = 0, to = 55, 5)) +
  ylab( "Weighted individuals" ) +
  xlab ("Generations") +
  scale_shape_manual(values = c(22,15, 17, 2, 19 )) +
  #scale_color_manual(values = c(c, dw, up))  +
  labs(shape = "lines", color ="sex") +
  panel_border() +
  theme(legend.position = c(0.8, 0.2),
        legend.direction = "horizontal",
        strip.background = element_rect(fill = "transparent") )

save_plot(filename = "demography_by_generation.png", plot = demo.plot, 
          base_aspect_ratio = 0.9, base_height = 8)


Demog %>% ggplot(aes(x=GEN, y=Census, group = interaction(selection, line), color = selection, shape =original.line) ) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ SEX, nrow = 2) +
  scale_x_discrete(breaks = seq(from = 0, to = 55, 5), labels = seq(from = 0, to = 55, 5)) +
  theme_bw() +
  xlab ("Generations") 


colors = c(c ="#CC79A7", dw = "#D55E00", up = "#0072B2")


mean.plot <- MeanP49 %>% ggplot(aes(x=GEN, y=P49, group = interaction(selection, line), color = selection, shape =original.line) ) +
  geom_point(size =3) +
  geom_line() +
  geom_point(data = hand.plot, shape = 16, color = "darkgrey", size = 5, alpha = 0.8) +
  facet_wrap(~ SEX, nrow = 2) +
  scale_x_discrete(breaks = seq(from = 0, to = 55, 5), labels = seq(from = 0, to = 55, 5)) +
  ylab( "Weight at 49 days (g)" ) +
  xlab ("Generations") +
  scale_shape_manual(values = c(22,15, 17, 2, 19 )) +
  scale_color_manual(values = c(c, dw, up))  +
  labs(shape = "lines") +
  panel_border() +
  theme(legend.position = "bottom", 
        strip.background = element_rect(fill = "transparent"))  
save_plot(filename = "mean_P49_by_generation.png", plot = mean.plot, 
          base_aspect_ratio = 0.9, base_height = 8)

