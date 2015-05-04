library(plyr)
library(dplyr)
library(magrittr)
library(lme4)
library(ggplot2)
library(tidyr)
library(evolqg)
library(readr)

raw.data <- tbl_df(read_csv("~/Dropbox/Ratones (1)//RRatones/Análises/Ratabase_Main.csv"))
raw.data %<>% mutate(treatment = LIN, strain = LIN)

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
raw.data %<>% filter(P49 < 50)

ggplot(raw.data, aes(P49, group = strain, color = interaction(treatment, strain))) + geom_histogram() + facet_grid(strain~treatment)
ggplot(filter(raw.data, strain == 'control'), aes(SEX, P49, color= SEX)) + geom_violin() + geom_jitter()

raw.main.data <- dlply(raw.data, .(treatment, strain), tbl_df)

current.data <- raw.main.data[[1]] 
makeMainData <- function (current.data) {
  x = vector("list", 8)
  x[[1]] <- select(current.data, c(ID:TAKE, strain, treatment))
  x[[2]] <- select(current.data, IS_PM:BA_OPI)
  x[[3]] <- unique(select(current.data, c(ID:P49, strain, treatment)))
  x[[4]] <- ddply(select(current.data, c(ID, IS_PM:BA_OPI)), .(ID), numcolwise(mean))[,-1]
  rownames(x[[4]]) <- x[[3]]$ID
  names(x)[1:4] <- c('info.raw', 'ed.raw', 'info', 'ed')
  x[[5]] <- CalcRepeatability(current.data$ID, ind.data = x$ed.raw)
  names(x)[5] <- 'reps'
  x[[6]] <- lm(as.matrix(x$ed) ~ x$info$SEX)
  x[[7]] <- CalculateMatrix(x[[6]])
  x[[8]] <- colMeans(x$ed)
  names(x)[6:8] <- c('model', 'cov.matrix', 'ed.means')
  return(x)
}
main.data <- llply(raw.main.data, makeMainData)