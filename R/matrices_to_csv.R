trait.names = colnames(raw.data)[10:44]
g_matrices =llply(g_models, function(x) alply(x$Gs, 1) ) 
p_matrices =llply(g_models, function(x) alply(x$Ps, 1) ) 

matrices.to.csv = function(x, line.name, type)  
{
  y = x 
  for (i in 1:length(x)) 
  {
    y[[i]] = as.data.frame(mutate(tbl_df(x[[i]]), matrix = paste(line.name, type, i, sep = "_"), trait = trait.names) )
  } 
  z = tbl_df(ldply(y))
  z = z %>% select(matrix, trait, V1:V35)
  colnames(z) = c("Matrix", "", trait.names)
  return(z)
}

############################################################################################
#                               G-Matrices                                                 #
############################################################################################

All.t.G = matrices.to.csv(p_matrices$control.t, line.name = "control_t", type = "G" )
write.csv(All.t.G, file = "data/SI_Matrices/G/All_control_t_posterior_G.csv")

All.s.G = matrices.to.csv(p_matrices$downwards.s, line.name = "downwards_s", type = "G" )
write.csv(All.s.G, file = "data/SI_Matrices/G/All_downwards_s_posterior_G.csv")

All.h.G = matrices.to.csv(p_matrices$downwards.h, line.name = "downwards_h", type = "G" )
write.csv(All.h.G, file = "data/SI_Matrices/G/All_downwards_h_posterior_G.csv")

All.sp.G = matrices.to.csv(p_matrices$`upwards.s'`, line.name = "upwards_s'", type = "G" )
write.csv(All.sp.G, file = "data/SI_Matrices/G/All_upwards_s'_posterior_G.csv")

All.hp.G = matrices.to.csv(p_matrices$`upwards.h'`, line.name = "upwards_h'", type = "G" )
write.csv(All.hp.G, file = "data/SI_Matrices/G/All_upwards_h'_posterior_G.csv")

############################################################################################
#                               P-Matrices                                                 #
############################################################################################

All.t.P = matrices.to.csv(p_matrices$control.t, line.name = "control_t", type = "P" )
write.csv(All.t.P, file = "data/SI_Matrices/P/All_control_t_posterior_P.csv")

All.s.P = matrices.to.csv(p_matrices$downwards.s, line.name = "downwards_s", type = "P" )
write.csv(All.s.P, file = "data/SI_Matrices/P/All_downwards_s_posterior_P.csv")

All.h.P = matrices.to.csv(p_matrices$downwards.h, line.name = "downwards_h", type = "P" )
write.csv(All.h.P, file = "data/SI_Matrices/P/All_downwards_h_posterior_P.csv")

All.sp.P = matrices.to.csv(p_matrices$`upwards.s'`, line.name = "upwards_s'", type = "P" )
write.csv(All.sp.P, file = "data/SI_Matrices/P/All_upwards_s'_posterior_P.csv")

All.hp.P = matrices.to.csv(p_matrices$`upwards.h'`, line.name = "upwards_h'", type = "P" )
write.csv(All.hp.P, file = "data/SI_Matrices/P/All_upwards_h'_posterior_P.csv")

matrices.to.csv2 = function(x, line.name, type)  
{
  z = ldply(x) %>% tbl_df() %>% mutate(., matrix = paste(line.name, type, sep = "_"), trait = rep(trait.names, length(x)), statistic = .id) 
  
  z %<>% select(., matrix, statistic, trait, 2:36)
  colnames(z) = c("Matrix", "Statistic" ,"", trait.names)
  return(z)
}
 
############################################################################################
#                               G-Matrices                                                 #
############################################################################################

Stats.t.G = matrices.to.csv2(ic.mean.G.matrices$control.t, line.name = "control_t", type = "G" )
write.csv(Stats.t.G, file = "data/SI_Matrices/G/Stats_control_t_distribuition_G.csv")

Stats.s.G = matrices.to.csv2(ic.mean.G.matrices$downwards.s, line.name = "downwards_s", type = "G" )
write.csv(Stats.s.G, file = "data/SI_Matrices/G/Stats_downwards_s_distribuition_G.csv")

Stats.h.G = matrices.to.csv2(ic.mean.G.matrices$downwards.h, line.name = "downwards_h", type = "G" )
write.csv(Stats.h.G, file = "data/SI_Matrices/G/Stats_downwards_h_distribuition_G.csv")

Stats.sp.G = matrices.to.csv2(ic.mean.G.matrices$`upwards.s'`, line.name = "upwards_s'", type = "G" )
write.csv(Stats.sp.G, file = "data/SI_Matrices/G/Stats_upwards_s'_distribuition_G.csv")

Stats.hp.G = matrices.to.csv2(ic.mean.G.matrices$`upwards.h'`, line.name = "upwards_h'", type = "G" )
write.csv(Stats.hp.G, file = "data/SI_Matrices/G/Stats_upwards_h'_distribuition_G.csv")
