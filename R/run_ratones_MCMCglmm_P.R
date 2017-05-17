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

x = main.data[[1]]
runMCMCmodelsRatones <- function (x) {
  x$full$P49 %<>% log
  trait_names = names(x$full)[c(8, 12:46)]
  full_data_scaled <- x$full
  scaled_traits = scale(x$full[,trait_names])
  vars = attributes(scaled_traits)$`scaled:scale`
  full_data_scaled[,trait_names] = scaled_traits
  
  value = paste("cbind(",
                paste(trait_names, collapse = ', '),
                ")", sep = '')

  fixed_effects = "trait:SEX + trait:AGE - 1"

  mcmc_formula = as.formula(paste(value, fixed_effects, sep = " ~ "))

  prior <- list(R = list(V = diag(num.traits), n = 2.5))
  ratones_model_corr <- MCMCglmm(mcmc_formula,
                            data = full_data_scaled,
                            rcov = ~us(trait):units,
                            family = rep("gaussian", num.traits),
                            nitt=103000, thin = 100, burnin = 3000,
                            prior = prior,
                            verbose = TRUE)
  corrPs = array(ratones_model_corr$VCV, dim = c(1000, num.traits, num.traits))
  corrP = apply(corrPs, 2:3, mean)

  varPs = outer(sqrt(vars), sqrt(vars))
  Ps = aaply(corrPs, 1, function(x) x * varPs)
  P = apply(Ps, 2:3, median)
  return(list(Ps = Ps,
              corrPs = corrPs,
              P = P,
              corrP = corrP,
              mcor = ratones_model_corr))
}

r_models = llply(main.data, runMCMCmodelsRatones, .parallel = TRUE)
new_names = names(main.data)
for(i in 1:length(r_models)) r_models[[i]]$line <- new_names[[i]]
#save(r_models, file = "Rdatas/ratonesMCMCmodels.RData")
load("Rdatas/ratonesMCMCmodels.RData")


