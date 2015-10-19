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


vectorCor <- function(x, y) t(Normalize(x)) %*% Normalize(y)

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

raw.data <- tbl_df(read_csv("./data/Ratabase_Main.csv"))
raw.data %<>% mutate(treatment = LIN, line = LIN)

#line colors
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
raw.data$line    <- gsub('t', 'control', raw.data$line)
raw.data$line    <- gsub('hp', 'h', raw.data$line)
raw.data$line    <- gsub('sp', 's', raw.data$line)

#Remove fat fucks
raw.data %<>% filter(P49 < 50) %>%
              filter(ID != 270160) %>% #outlier in biplot prcomp - control
              filter(ID != 270493) %>% #outlier in biplot prcomp - increase.S
              filter(ID != 270200) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270556) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270190)     #outlier in biplot prcomp - reduce.s


ggplot(raw.data, aes(P49, group = line, color = interaction(treatment, line))) + geom_histogram() + facet_grid(line~treatment)
ggplot(filter(raw.data, line == 'control'), aes(SEX, P49, color= SEX)) + geom_violin() + geom_jitter()

raw.main.data <- dlply(raw.data, .(treatment, line), tbl_df)

current.data <- raw.main.data[[4]]



makeMainData <- function (current.data) {
  x = vector("list", 11)
  current.data$AGE[is.na(current.data$AGE)] <- mean(current.data$AGE, na.rm = TRUE)
  x[[1]] <- select(current.data, c(ID:TAKE, line, treatment))
  x[[2]] <- select(current.data, c(P49, IS_PM:BA_OPI))
  x[[3]] <- unique(select(current.data, c(ID:P49, line, treatment)))
  x[[4]] <- ddply(select(current.data, c(ID, P49, IS_PM:BA_OPI)), .(ID), numcolwise(mean))
  set_row <- function(x) {rownames(x) <- x$ID; x[,-1]}
  x[[4]] <- set_row(x[[4]])
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed')
  x[[5]] <- CalcRepeatability(current.data$ID, ind.data = x$ed.raw[,-1])
  names(x)[5] <- 'reps'
  sex_age_lm <- lm(as.matrix(x$ed) ~ x$info$SEX + x$info$AGE)
  sex_age_res <- residuals(sex_age_lm)
  p49_traits_pls <- plsreg1(sex_age_res[,2:36], sex_age_res[,1])
  #x[[10]] <- Normalize(p49_traits_pls$reg.coefs[-1])
  x[[10]] <- p49_traits_pls$reg.coefs[-1]
  x$ed$P49 %<>% log
  x[[6]] <- lm(as.matrix(x$ed) ~ x$info$SEX)
  x[[7]] <- CalculateMatrix(x[[6]])
  x[[8]] <- colMeans(x$ed)
  x[[9]] <- tbl_df(cbind(x$info, x$ed))
  x[[11]] <- mean((apply(x[['ed']], 1, gm_mean)))
  names(x)[6:11] <- c('model', 'cov.matrix', 'ed.means', 'full', 'plsr', 'gm_mean')
  return(x)
}
main.data <- llply(raw.main.data, makeMainData)

full_data = ldply(main.data, function(x) x$full)
Wmat <- CalculateMatrix(lm(as.matrix(select(full_data, IS_PM:BA_OPI)) ~ full_data$SEX*full_data$LIN))

main.data %>% laply(function(x) x$plsr) %>% {. %*% t(.)}

m_full_data = melt(full_data, id.vars = names(full_data)[c(1:8, 10, 11)])

full_trait_plots = ggplot(m_full_data, aes(line, value, group = interaction(treatment, line), fill = treatment)) + geom_boxplot() + facet_wrap(~variable, scale = "free") + background_grid(major = 'y', minor = "none")
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/full_trait_plot.pdf", full_trait_plots, ncol = 6, nrow = 6)
