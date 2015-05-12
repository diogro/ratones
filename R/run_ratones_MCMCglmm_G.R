source('./R/read_ratones.R')
if(!require (MasterBayes)) {install.packages("MasterBayes"); library(MasterBayes)}

ped <- unique(read.csv("./data/pedigree.csv"))
#fixing pedigree so ID are not repeated
count_ID <- table(ped[,1])
ped <- ped[!(ped[,1] %in% names(which(count_ID > 1)) & is.na(ped[,2])),]
ped <- orderPed (ped)

num.traits = 36

x <- main.data[[1]]
full_data <- cbind(x$info[,-8], scale(x$ed))
full_data$animal <- full_data$ID

value = paste("cbind(",
              paste(names(full_data)[10:45], collapse = ', '),
              ")", sep = '')

fixed_effects = "trait:SEX - 1"

mcmc_formula = as.formula(paste(value, fixed_effects, sep = " ~ "))

prior <- list(R = list(V = diag(num.traits), 
                       n = 0.002),
              G = list(G1 = list(V = diag(num.traits) * 0.02, 
                                 n = num.traits + 1)))
ratones_model <- MCMCglmm(mcmc_formula, 
                          random = ~idh(trait):animal,
                          data = full_data,
                          rcov = ~idh(trait):units,
                          family = rep("gaussian", num.traits),
                          pedigree = ped,
                          nitt=200000, thin = 100, burnin = 100000,
                          prior = prior,
                          verbose = TRUE)

posterior.mode(ratones_model$VCV)[1:36]
diag(cov((full_data)[10:45]))
summary(ratones_model)
