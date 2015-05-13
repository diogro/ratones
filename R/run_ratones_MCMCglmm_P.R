source('./R/read_ratones.R')

num.traits = 36

x <- main.data[[1]]
x$ed.raw$P49 %<>% log
full_data <- cbind(x$info[,-8], scale(x$ed))

value = paste("cbind(",
              paste(names(select(full_data, P49:BA_OPI)), collapse = ', '),
              ")", sep = '')

fixed_effects = "trait:SEX - 1"

mcmc_formula = as.formula(paste(value, fixed_effects, sep = " ~ "))
 
# prior <- list(R = list(V = diag(num.traits), 
#                        n = 0.002),
#               G = list(G1 = list(V = diag(num.traits) * 0.1, 
#                                  n = num.traits + 1)))

prior1 <- list(R = list(V = diag(num.traits), n = 0.1))
ratones_model_1 <- MCMCglmm(mcmc_formula, 
#                           random = ~idh(trait):ID,
                          data = full_data,
                          rcov = ~us(trait):units,
                          family = rep("gaussian", num.traits),
                          nitt=100000, thin = 100, burnin = 10000,
                          prior = prior1,
                          verbose = TRUE)
prior2 <- list(R = list(V = diag(num.traits), n = num.traits + 1))
ratones_model_2 <- MCMCglmm(mcmc_formula, 
                          #                           random = ~idh(trait):ID,
                          data = full_data,
                          rcov = ~us(trait):units,
                          family = rep("gaussian", num.traits),
                          nitt=100000, thin = 100, burnin = 10000,
                          prior = prior2,
                          verbose = TRUE)

dim(ratones_model_1$Sol)

Ps_1 = array(ratones_model_1$VCV, dim = c(900, num.traits, num.traits))
Ps_1 = aaply(Ps_1, 1, cov2cor)
P1 = apply(Ps_1, 2:3, mean)

Ps_2 = array(ratones_model_2$VCV, dim = c(900, num.traits, num.traits))
Ps_2 = aaply(Ps_2, 1, cov2cor)
P2 = apply(Ps_2, 2:3, mean)

CalcR2(P1)
CalcR2(P2)
CalcR2(x$cov.matrix)
MatrixCompare(P1, P2)
MatrixCompare(P1, cov2cor(x$cov.matrix))
MatrixCompare(P2, cov2cor(x$cov.matrix))
ggplot(melt(data.frame(montecarlo = MonteCarloR2(x$cov.matrix, 36, iterations = 900),
                       weak = apply(Ps_1, 1, CalcR2), 
                       strong = apply(Ps_2, 1, CalcR2))), 
       aes(variable, value)) + geom_boxplot()
