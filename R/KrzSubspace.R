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
observed = HPDinterval(as.mcmc(MCMC.H.val), prob = 0.95)


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
null = HPDinterval(as.mcmc(null_MCMC.H.val), prob = 0.95)
krz_subspace_plot = rbind(cbind(rank = 1:35, as.data.frame(observed), type = "Observed"), 
      cbind(rank = 1:35, as.data.frame(null), type = "Randomised")) %>%
  mutate(mean = (upper + lower) / 2) %>%
  ggplot(aes(x = rank, y = mean, linetype = type, color = type)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
  labs(y = "Eigenvalues of H", x = "Eigenvectors of H") + background_grid(major = 'x', minor = "none") + 
  panel_border() + theme(panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.2))+ 
  theme(legend.position = c(0.75, 0.9), text = element_text(size = 18)) + 
  scale_colour_grey("", start = 0, end = 0.6) + scale_linetype(guide = "none")
  krz_subspace_plot
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/krz_subspace_(figure2?).png", krz_subspace_plot, base_aspect_ratio = 1.3, base_height = 4.8)

dimnames(p_matrices)[[3]] = Gnames

PlotRSprojection_rata <- function( rs_proj = rs_projection, cov.matrix.array = p_matrices, num_pc = 8, p = 0.95, ncols = 5, label = "Phenotypic Variance")
{
n <- dim(cov.matrix.array)[[1]]
evecs = t(rs_proj$eig.R$vectors)
R_eigen_projection <- apply(cov.matrix.array, 3:4, function(mat) diag(evecs %*% 
                                                                        mat %*% t(evecs)))
R_eigen_projection <- aperm(R_eigen_projection, c(3, 2, 1))
colnames(R_eigen_projection) <- dimnames(cov.matrix.array)[[3]]
HPD.int <- aaply(R_eigen_projection, 3, function(proj) HPDinterval(as.mcmc(proj), 
                                                                   prob = p))
HPD.int <- aperm(HPD.int, c(2, 3, 1))
dimnames(HPD.int) = list(dimnames(cov.matrix.array)[[3]], 
                         NULL, NULL)
dat = adply(HPD.int, 1:3)
names(dat) = c("Population", "interval", "trait", "value")
dat$trait <- as.numeric(dat$trait)
if(!is.null(num_pc))  dat %<>% filter(trait <=num_pc)
dat$trait <- as.factor(dat$trait)
dat = dcast(dat, as.formula("Population+trait~interval"))
names(dat) = c("Population", "trait", "lower", "upper")
dat$trait = paste0("E", dat$trait)
order.list = paste0("E", 1:num_pc)
dat$trait = factor(dat$trait, levels = order.list)
dat$Population <- as.factor(dat$Population)
dat$Population = factor(dat$Population, levels = c("control.t", "upwards.h'", "upwards.s'","downwards.h", "downwards.s"))
dat$mean = rowMeans(cbind(dat$upper, dat$lower))
myPalette <- viridis(5)

plot = ggplot(dat, aes_string( colour= "Population", x = "Population", y = "mean")) + geom_point() + 
  geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), size = 1 ) + 
  scale_color_manual(values = myPalette, name = "Line", labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) +
  theme(axis.text.x = element_text(hjust = 1)) + background_grid() + 
  theme(legend.position = c(0.85, 0.1), text = element_text(size = 18)) +
  scale_x_discrete("Lines", labels = c("control.t" = "t","upwards.h'" = "h'", "upwards.s'" = "s'","downwards.h" = "h","downwards.s" = "s")) +
  facet_wrap(~trait, ncol = ncols, scales = "free_y") + panel_border() + ylab(label) 

return(plot)
}

rs_projection = RSProjection(p_matrices)
rs_projection_plot_full = PlotRSprojection_rata(rs_proj = rs_projection, cov.matrix.array = p_matrices, num_pc = 35, p = 0.95, ncols = 5)
rs_projection_plot_full = rs_projection_plot_full +  theme(legend.position = "bottom")
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/rs_projection_(full_SI).png", rs_projection_plot_full, base_aspect_ratio = 1.1, base_height = 1.8, ncol = 5, nrow = 7)

rs_projection_plot = PlotRSprojection_rata(rs_proj = rs_projection, p_matrices, p = 0.95, num_pc = 8, ncols = 3)
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/rs_projection_(figure2?).png", rs_projection_plot, base_aspect_ratio = 1.3, base_height = 4.8)

figure_2_review = plot_grid(krz_subspace_plot, rs_projection_plot, labels = c("A", "B"), ncol = 2, rel_widths = c(1,1.3))
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/figure2_review.png", figure_2_review, base_aspect_ratio = 1.3, base_height = 4.8, ncol = 2)

#G matrices

g_matrices = laply(g_models, function(x) x$Gs)
g_matrices = aperm(g_matrices, c(3, 4, 1, 2))
dimnames(g_matrices)[[3]] <- dimnames(p_matrices)[[3]]
rs_projection_plot_g = PlotRSprojection_rata(rs_proj = rs_projection, g_matrices, num_pc = 8, p = 0.95, ncols = 3, label = "Genetic Variance")
save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/rs_projection_g.png", rs_projection_plot_g, base_aspect_ratio = 1.3, base_height = 4.8)



