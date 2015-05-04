library(plyr)
library(dplyr)
library(magrittr)
library(lme4)
library(ggplot2)
library(tidyr)
library(evolqg)
library(readr)

raw.data <- tbl_df(read_csv("~/Dropbox/Ratones (1)//RRatones/Análises/Ratabase_Main.csv"))
raw.data %<>% mutate(treatment = LIN,
                     strain = LIN)

# Change weird labels
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

ggplot(raw.data, aes(P49, group = strain, color = interaction(treatment, strain))) + geom_histogram() + facet_grid(treatment~strain)
ggplot(filter(raw.data, strain == 'control'), aes(SEX, P49, color= SEX)) + geom_violin() + geom_jitter()
