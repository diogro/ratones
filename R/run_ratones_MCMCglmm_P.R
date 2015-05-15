source('./R/read_ratones.R')
if(!require(doMC)) {install.packages('doMC'); library(doMC)}

registerDoMC(5)
num.traits = 36

find_CI = function(x, prob = 0.95){
  n = length(x)
  xs = sort(x)
  nint = floor(prob*n)
  lowest_int = abs(xs[n] - xs[1])
  for(i in 1:(n-nint)){
    current_int = abs(xs[i] - xs[i+nint])
    if(current_int <= lowest_int){
      lowest_int = current_int
      pos = i
    }
  }
  return(c(xs[pos], xs[pos+nint]))
}

find_CI_upper = function(x, prob = 0.95){
  n = length(x)
  xs = sort(x)
  nint = floor(prob*n)
  lowest_int = abs(xs[n] - xs[1])
  for(i in 1:(n-nint)){
    current_int = abs(xs[i] - xs[i+nint])
    if(current_int <= lowest_int){
      lowest_int = current_int
      pos = i
    }
  }
  return(c(upper = xs[pos+nint]))
}

find_CI_lower = function(x, prob = 0.95){
  n = length(x)
  xs = sort(x)
  nint = floor(prob*n)
  lowest_int = abs(xs[n] - xs[1])
  for(i in 1:(n-nint)){
    current_int = abs(xs[i] - xs[i+nint])
    if(current_int <= lowest_int){
      lowest_int = current_int
      pos = i
    }
  }
  return(c(lower = xs[pos]))
}
runMCMCmodelsRatones <- function (x) {
  x$ed.raw$P49 %<>% log
  full_data_scaled <- cbind(x$info[,-8], scale(x$ed))
  
  value = paste("cbind(",
                paste(names(select(full_data_scaled, P49:BA_OPI)), collapse = ', '),
                ")", sep = '')
  
  fixed_effects = "trait:SEX + trait:AGE - 1"
  
  mcmc_formula = as.formula(paste(value, fixed_effects, sep = " ~ "))
  
  prior <- list(R = list(V = diag(num.traits), n = 2))
  ratones_model_corr <- MCMCglmm(mcmc_formula, 
                            data = full_data_scaled,
                            rcov = ~us(trait):units,
                            family = rep("gaussian", num.traits),
                            nitt=13000, thin = 100, burnin = 3000,
                            prior = prior,
                            verbose = TRUE)
  corrPs = array(ratones_model_corr$VCV, dim = c(100, num.traits, num.traits))
  corrPs = aaply(corrPs, 1, cov2cor)
  corrP = apply(corrPs, 2:3, mean)
  
  full_data <- cbind(x$info[,-8], x$ed)
  prior <- list(R = list(V = diag(num.traits), n = 0.01))
  ratones_model_var <- MCMCglmm(mcmc_formula, 
                            data = full_data,
                            rcov = ~idh(trait):units,
                            family = rep("gaussian", num.traits),
                            nitt=13000, thin = 100, burnin = 3000,
                            prior = prior,
                            verbose = TRUE)
  varPs = aaply(ratones_model_var$VCV, 1, function(x) outer(sqrt(x), sqrt(x)))
  Ps = corrPs * varPs
  P = apply(Ps, 2:3, mean)
  return(list(Ps = Ps,
              corrPs = corrPs, 
              P = P,
              corrP = corrP,
              mcor = ratones_model_corr,
              mvar = ratones_model_var))
}

#r_models = llply(main.data, runMCMCmodelsRatones, .parallel = TRUE)
#save(r_models, file = "Rdatas/ratonesMCMCmodels.RData")
load("Rdatas/ratonesMCMCmodels.RData")

