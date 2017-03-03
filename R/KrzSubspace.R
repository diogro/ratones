library(coda)
library(MCMCglmm)
library(MasterBayes)

p_matrices = laply(g_models, function(x) x$Ps)
p_matrices = aperm(p_matrices, c(3, 4, 1, 2))
Hs = alply(p_matrices, 4, function(x) alply(x, 3)) %>% llply(function(x) KrzSubspace(x)$H)
avgH = Reduce("+", Hs)/length(Hs)
avgH.vec <- eigen(avgH)$vectors
MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))

# confidence intervals for variation in shared subspace directions
observed = HPDinterval(as.mcmc(MCMC.H.val))


ped <- unique(read.csv("./data/pedigree.csv"))
#fixing pedigree so ID are not repeated
count_ID <- table(ped[,1])
ped <- ped[!(ped[,1] %in% names(which(count_ID > 1)) & is.na(ped[,2])),]
ped = orderPed(ped)
names(ped) = c("id", "dam", "sire")

n = 35
m = 5
MCMCsamp = 1000
traitnames = names(dplyr::select(full_data, IS_PM:BA_OPI))
Gnames = names(g_models)
rand.g_models <- array(NA ,c(n,n,m,MCMCsamp))
dimnames(rand.g_models) <- list(traitnames,traitnames,Gnames)
IDs = llply(main.data, function(x) x$info$ID)
n_IDs = laply(IDs, length)
filter_ped = function(x, i) x[rownames(x) %in% IDs[[i]],]
for (i in 1:MCMCsamp){
  ct.bv <- filter_ped(rbv(ped, g_models[[1]]$Ps[i,,]), 1)
  dh.bv <- filter_ped(rbv(ped, g_models[[2]]$Ps[i,,]), 2)
  ds.bv <- filter_ped(rbv(ped, g_models[[3]]$Ps[i,,]), 3)
  uh.bv <- filter_ped(rbv(ped, g_models[[4]]$Ps[i,,]), 4)
  us.bv <- filter_ped(rbv(ped, g_models[[5]]$Ps[i,,]), 5)
  a.pop <- cumsum(n_IDs)
  pop.bv <- rbind(ct.bv,dh.bv,ds.bv,uh.bv,us.bv)
  rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
  rand.g_models[,,1,i] <- cov(rand.pop.bv[1:a.pop[1],])
  rand.g_models[,,2,i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2],])
  rand.g_models[,,3,i] <- cov(rand.pop.bv[(a.pop[2] + 1):a.pop[3],])
  rand.g_models[,,4,i] <- cov(rand.pop.bv[(a.pop[3] + 1):a.pop[4],])
  rand.g_models[,,5,i] <- cov(rand.pop.bv[(a.pop[4] + 1):a.pop[5],])
}

null_Hs = alply(rand.g_models, 4, function(x) alply(x, 3)) %>% llply(function(x) KrzSubspace(x)$H)
null_avgH = Reduce("+", null_Hs)/length(null_Hs)
null_avgH.vec <- eigen(null_avgH)$vectors
null_MCMC.H.val = laply(null_Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
null = HPDinterval(as.mcmc(null_MCMC.H.val))
rbind(cbind(rank = 1:35, as.data.frame(observed), type = "observed"), 
      cbind(rank = 1:35, as.data.frame(null), type = "null")) %>%
  mutate(mean = (upper + lower) / 2) %>%
  ggplot(aes(x = rank, y = mean, linetype = type, shape = type)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5))

