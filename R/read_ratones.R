if(!require(plyr)) {install.packages('plyr'); library(plyr)}
if(!require(dplyr)) {install.packages('dplyr'); library(dplyr)}
if(!require(magrittr)) {install.packages('magrittr'); library(magrittr)}
if(!require(lme4)) {install.packages('lme4'); library(lme4)}
if(!require(ggplot2)) {install.packages('ggplot2'); library(ggplot2)}
if(!require(tidyr)) {install.packages('tidyr'); library(tidyr)}
if(!require(MCMCglmm)) {install.packages('MCMCglmm'); library(MCMCglmm)}
if(!require(reshape2)) {install.packages('reshape2'); library(reshape2)}
if(!require(evolqg)) {devtools::install_github('diogro/evolqg'); library(evolqg)}
if(!require(readr)) {devtools::install_github('hadley/readr'); library(readr)}
if(!require(cowplot)) {install.packages('cowplot'); library(cowplot)}
if(!require(plsdepot)) {install.packages('plsdepot'); library(plsdepot)}
if(!require(xtable)) {install.packages('xtable'); library(xtable)}

setwd("~/projects/ratones/")

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
raw.data %<>% mutate(selection = LIN, line = LIN)

# #line colors
c = "#CC79A7"
dw = "#D55E00"
up = "#0072B2"

#Line order
#lines = c("t", "h", "s", "h'", "s'")
lines = c("t", "h'", "s'", "h", "s")

#line colors
# c = rgb(255, 0, 255, maxColorValue = 255)
# h = rgb(0, 0, 128, maxColorValue = 255)
# s = rgb(0, 128, 128, maxColorValue = 255)

# Change weird labels
names(raw.data) <- gsub("-", "_", names(raw.data))
raw.data$selection <- gsub('t', 'control', raw.data$selection)
raw.data$selection <- gsub('hp', 'upwards', raw.data$selection)
raw.data$selection <- gsub('h', 'downwards', raw.data$selection)
raw.data$selection <- gsub('sp', 'upwards', raw.data$selection)
raw.data$selection <- gsub("\\bs\\b", 'downwards', raw.data$selection, perl = TRUE)
raw.data$original_line <- raw.data$line
raw.data$original_line <- gsub('hp', "h'", raw.data$original_line)
raw.data$original_line <- gsub('sp', "s'", raw.data$original_line)
raw.data$line    <- gsub('hp', "h'", raw.data$line)
raw.data$line    <- gsub('sp', "s'", raw.data$line)
raw.data$original_line <- factor(raw.data$original_line, levels = lines)
raw.data$line <- factor(raw.data$line, levels = lines)

#raw.data$line    <- gsub('t', 'control', raw.data$line)


#Remove fat fucks
raw.data %<>% filter(P49 < 50) %>%
              filter(ID != 270160) %>% #outlier in biplot prcomp - control
              filter(ID != 270493) %>% #outlier in biplot prcomp - increase.S
              filter(ID != 270200) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270556) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270190)     #outlier in biplot prcomp - reduce.s

raw.main.data <- dlply(raw.data, .(selection), tbl_df)

current.data <- raw.main.data[[2]]

makeMainData <- function (current.data) {
  x = vector("list", 11)
  current.data$AGE[is.na(current.data$AGE)] <- mean(current.data$AGE, na.rm = TRUE)
  x[[1]] <- select(current.data, c(ID:TAKE, line, original_line, selection))
  x[[2]] <- select(current.data, IS_PM:BA_OPI)
  x[[3]] <- unique(select(current.data, c(ID:P49, line, original_line, selection)))
  x[[4]] <- ddply(select(current.data, ID, IS_PM:BA_OPI), .(ID), numcolwise(mean))
  set_row <- function(x) {rownames(x) <- x$ID; x[,-1]}
  x[[4]] <- set_row(x[[4]])
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed')
  x[[5]] <- CalcRepeatability(current.data$ID, ind.data = x$ed.raw[,-1])
  names(x)[5] <- 'reps'
  if(length(unique(x$info$original_line))==1)   sex_age_lm <- lm(as.matrix(x$ed) ~ x$info$SEX + x$info$AGE)
  else sex_age_lm <- lm(as.matrix(x$ed) ~ x$info$original_line*x$info$SEX*x$info$AGE)
  sex_age_res <- residuals(sex_age_lm)
  p49_traits_pls <- plsreg1(sex_age_res[,2:35], sex_age_res[,1])
  #x[[10]] <- Normalize(p49_traits_pls$reg.coefs[-1])
  x[[10]] <- p49_traits_pls$reg.coefs[-1]
  #x$ed$P49 %<>% log
  if(length(unique(x$info$original_line))==1) x[[6]] <- lm(as.matrix(x$ed) ~ x$info$SEX + x$info$AGE)
  else x[[6]] <- lm(as.matrix(x$ed) ~ x$info$original_line*x$info$SEX*x$info$AGE)
  x[[7]] <- CalculateMatrix(x[[6]])
  x[[8]] <- colMeans(x$ed)
  x[[9]] <- tbl_df(cbind(x$info, x$ed))
  x[[11]] <- mean((apply(x[['ed']], 1, gm_mean)))
  names(x)[6:11] <- c('model', 'cov.matrix', 'ed.means', 'full', 'plsr', 'gm_mean')
  return(x)
}
main.data <- llply(raw.main.data, makeMainData)

full_data = ldply(main.data, function(x) x$full)
Wmat <- CalculateMatrix(lm(as.matrix(select(full_data, IS_PM:BA_OPI)) ~ full_data$AGE*full_data$SEX*full_data$LIN))

r_models = llply(main.data, function(x) BayesianCalculateMatrix(x$model, samples = 500, nu = 3))
for(i in 1:length(r_models)){
  r_models[[i]]$line <- names(r_models)[1]
}

for(i in 1:length(r_models)){
  if(any(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)){
    print("non symmetric")
    x = r_models[[i]]$Ps[(which(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)),,]
    r_models[[i]]$Ps[(which(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)),,] <- (x + t(x))/2
  }
}


res = residuals(main.data$downwards$model) %*% eigen(r_models$downwards$MAP)$vectors
cov(res)
res_l = data.frame(res, main.data$downwards$info)
ggplot(res_l, aes(X1, X2, color = line)) + geom_point()
