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
if(!require(viridis)) {install.packages('viridis'); library(viridis)}

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
raw.data %<>% filter(P49 < 50)  %>%
              filter(ID != 270160) %>% #outlier in biplot prcomp - control
              filter(ID != 270493) %>% #outlier in biplot prcomp - increase.S
              filter(ID != 270200) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270556) %>% #outlier in biplot prcomp - reduce.h
              filter(ID != 270190)     #outlier in biplot prcomp - reduce.s

raw.main.data <- dlply(raw.data, .(selection, line), tbl_df)

current.data <- raw.main.data[[2]]

makeMainData <- function (current.data) {
  x = vector("list", 9)
  names(x) <- c('info.raw', 'ed.raw', 'info', 'ed', 'reps', 'model', 'ed.means', 'full', 'gm_mean')
  current.data$AGE[is.na(current.data$AGE)] <- mean(current.data$AGE, na.rm = TRUE)
  x$info.raw = dplyr::select(current.data, c(ID:TAKE, line, original_line, selection))
  x$ed.raw   = dplyr::select(current.data, ID, IS_PM:BA_OPI)
  x$info     = arrange(unique(dplyr::select(current.data, c(ID:P49, line, original_line, selection))), ID)
  x$ed       = arrange(ddply(dplyr::select(current.data, ID, IS_PM:BA_OPI), .(ID), numcolwise(mean)), ID)
  x$full <- tbl_df(inner_join(x$info, x$ed, by = "ID"))
  set_row <- function(x) {rownames(x) <- x$ID; x[,-1]}
  x$ed <- set_row(x$ed)
  x$reps <- CalcRepeatability(x$ed.raw$ID, ind.data = x$ed.raw[,-1])
  if(length(unique(x$info$original_line))==1) x$model <- lm(as.matrix(x$ed) ~ x$info$SEX + x$info$AGE)
  else x$model <- lm(as.matrix(x$ed) ~ x$info$original_line + x$info$SEX + x$info$AGE)
  x$ed.means <- colMeans(x$ed)
  x$gm_mean <- mean((apply(x$'ed', 1, gm_mean)))
  return(x)
}
main.data <- llply(raw.main.data, makeMainData)

full_data = ldply(main.data, function(x) x$full)
Wmat <- CalculateMatrix(lm(as.matrix(dplyr::select(full_data, IS_PM:BA_OPI)) ~ full_data$AGE*full_data$SEX*full_data$LIN))

r_models = llply(main.data, function(x) BayesianCalculateMatrix(x$model, samples = 500, nu = 3))
for(i in 1:length(r_models)){
  r_models[[i]]$line <- names(r_models)[i]
}

for(i in 1:length(r_models)){
  if(any(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)){
    print("non symmetric")
    x = r_models[[i]]$Ps[(which(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)),,]
    r_models[[i]]$Ps[(which(aaply(r_models[[i]]$Ps, 1, isSymmetric)==FALSE)),,] <- (x + t(x))/2
  }
}

folder = "t"
readMatLab <- function(folder){
    x = readMat(paste0("./Gmatlab/", folder, "/Posterior_mean.mat"))
    G = x$posterior.mean[,,1]$G
    P = x$posterior.mean[,,1]$P
    Gs = aperm(x$posterior.mean[,,1]$Gs, c(3, 1, 2))
    Ps = aperm(x$posterior.mean[,,1]$Ps, c(3, 1, 2))
    list(P = P, Ps = Ps, G = G, Gs = Gs)
    }
g_models = list(control.t = readMatLab("t"),
     downwards.h = readMatLab("h"),
     downwards.s = readMatLab("s"),
     "upwards.h'" = readMatLab("hp"),
     "upwards.s'" = readMatLab("sp"))
for(i in 1:length(g_models)){ g_models[[i]]$line <- names(g_models)[i] }

x = readMat("./Gmatlab/ratones/Posterior_mean.mat")
G = x$posterior.mean[,,1]$G
P = x$posterior.mean[,,1]$P

delta_Z = full_data %>% filter(selection == "upwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans -
    full_data %>% filter(selection == "downwards") %>% dplyr::select(IS_PM:BA_OPI) %>% colMeans
