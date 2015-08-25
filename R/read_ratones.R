if(!require(plyr)) {install.packages('plyr'); library(plyr)}
if(!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if(!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
if(!require(lme4)) {install.packages('lme4'); library(lme4)}
if(!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if(!require(tidyr)) {install.packages('tidyr'); library(tidyr)}
if(!require(MCMCglmm)) {install.packages('MCMCglmm'); library(MCMCglmm)}
if(!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
if(!require(evolqg)) {devtools::install_github('lem-usp/evolqg'); library(evolqg)}
if(!require(readr)) {devtools::install_github('hadley/readr'); library(readr)}
if(!require(cowplot)) {install.packages('cowplot'); library(cowplot)}
if(!require(plsdepot)) {install.packages('plsdepot'); library(plsdepot)}


raw.data <- tbl_df(read_csv("./data/Ratabase_Main.csv"))
raw.data %<>% mutate(treatment = LIN, strain = LIN)

#Strain colors
c = "#CC79A7"
h = "#D55E00"
s = "#0072B2"

# Change weird labels
names(raw.data) <- gsub("-", "_", names(raw.data))
raw.data$treatment <- gsub('t', 'control', raw.data$treatment)
raw.data$treatment <- gsub('hp', 'increase', raw.data$treatment)
raw.data$treatment <- gsub('h', 'reduce', raw.data$treatment)
raw.data$treatment <- gsub('sp', 'increase', raw.data$treatment)
raw.data$treatment <- gsub("\\bs\\b", 'reduce', raw.data$treatment, perl = TRUE)
raw.data$strain    <- gsub('t', 'control', raw.data$strain)
raw.data$strain    <- gsub('hp', 'h', raw.data$strain)
raw.data$strain    <- gsub('sp', 's', raw.data$strain)

#Remove fat fucks
raw.data %<>% filter(P49 < 50) %>%
              filter(ID != 270160) %>% #outlier in biplot prcomp - control
              filter(ID != 270493) %>% #outlier in biplot prcomp - increase.S
              filter(ID != 270200) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270556) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270190)     #outlier in biplot prcomp - reduce.s


ggplot(raw.data, aes(P49, group = strain, color = interaction(treatment, strain))) + geom_histogram() + facet_grid(strain~treatment)
ggplot(filter(raw.data, strain == 'control'), aes(SEX, P49, color= SEX)) + geom_violin() + geom_jitter()

raw.main.data <- dlply(raw.data, .(treatment, strain), tbl_df)

current.data <- raw.main.data[[4]]



makeMainData <- function (current.data) {
  x = vector("list", 10)
  current.data$AGE[is.na(current.data$AGE)] <- mean(current.data$AGE, na.rm = TRUE)
  x[[1]] <- select(current.data, c(ID:TAKE, strain, treatment))
  x[[2]] <- select(current.data, c(P49, IS_PM:BA_OPI))
  x[[3]] <- unique(select(current.data, c(ID:P49, strain, treatment)))
  x[[4]] <- ddply(select(current.data, c(ID, P49, IS_PM:BA_OPI)), .(ID), numcolwise(mean))
  set_row <- function(x) {rownames(x) <- x$ID; x[,-1]}
  x[[4]] <- set_row(x[[4]])
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed')
  x[[5]] <- CalcRepeatability(current.data$ID, ind.data = x$ed.raw[,-1])
  names(x)[5] <- 'reps'
  x$ed$P49 %<>% log
  x[[6]] <- lm(as.matrix(x$ed) ~ x$info$SEX)
  x[[7]] <- CalculateMatrix(x[[6]])
  x[[8]] <- colMeans(x$ed)
  x[[9]] <- tbl_df(cbind(x$info, x$ed))
  sex_age_lm <- lm(as.matrix(x$ed) ~ x$info$SEX + x$info$AGE)
  sex_age_res <- residuals(sex_age_lm)
  p49_traits_pls <- plsreg1(sex_age_res[,2:36], sex_age_res[,1])
  x[[10]] <- Normalize(p49_traits_pls$reg.coefs[-1])
  names(x)[6:10] <- c('model', 'cov.matrix', 'ed.means', 'full', 'plsr')
  return(x)
}
main.data <- llply(raw.main.data, makeMainData)

full_data = ldply(main.data, function(x) x$full)
Wmat <- CalculateMatrix(lm(as.matrix(select(full_data, IS_PM:BA_OPI)) ~ full_data$SEX*full_data$LIN))

main.data %>% laply(function(x) x$plsr) %>% {. %*% t(.)}

m_full_data = melt(full_data, id.vars = names(full_data)[c(1:8, 10, 11)])
full_trait_plots = ggplot(m_full_data, aes(strain, value, group = interaction(treatment, strain), fill = treatment)) + geom_boxplot() + facet_wrap(~variable, scale = "free")
ggsave("~/Desktop/full_trait_plot.pdf", full_trait_plots)
