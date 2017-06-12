#' Run the Krzanowski subspace comparison
#'
#' Compares the H matrix eigenvalues between observed and randomized datasets,
#' using simulated beeding values from the pedigree.
#' @param matrix_array array of posterior matrices with dimensions traits x traits x populations x MCMCsamples
#' @param models list of models, like ratones_models. Ignored if matrix_array is present.
#' @param type use P or G matrices from the models list. Ignored if matrix_array is present.
#' @param ped pedigree for null model simulation
#' @param IDs list of IDs for the individuals in each population. Each element must correspond to a populations, and IDs must be in the pedigree.
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the posterior intervals for the eigenvalues of H.
#' @return A data.frame with the observed and randomized intervals for the eigenvalues of H
#' @examples
#' \dontrun{
#' out = runKrzSubspace(models = ratones_models, type = "P", ped = ratones_ped$ped)
#' #figure 2, Panel A:
#' krz_subspace_plot = ggplot(out, aes(x = rank, y = mean, linetype = type, color = type)) +
#'   geom_point(position = position_dodge(width = 0.5)) +
#'   geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
#'   labs(y = "Eigenvalues of H", x = "Eigenvectors of H") + background_grid(major = 'x', minor = "none") +
#'   panel_border() + theme(panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size =0.2))+
#'   theme(legend.position = c(0.75, 0.9), text = element_text(size = 18)) +
#'   scale_colour_grey("", start = 0, end = 0.6) + scale_linetype(guide = "none")
#' #Panel B can be created in evolqg (but will take a long time!!):
#' matrix_array = laply(ratones_models, function(x) x$Ps)
#' matrix_array = aperm(matrix_array, c(3, 4, 1, 2))
#' dimnames(matrix_array)[[3]] = names(ratones_models)
#' rs_projection = RSProjection(matrix_array)
#' PlotRSprojection(rs_proj = rs_projection, cov.matrix.array = matrix_array)
#'}
runKrzSubspace = function(matrix_array = NULL, models = NULL, type = c("P", "G"), ped = NULL,
                          IDs = llply(ratones, function(x) x$info$ID), prob = 0.95){
  type = match.arg(type)
  if(is.null(matrix_array)){
    if(type == "G") matrix_array = laply(models, function(x) x$Gs)
    else matrix_array = laply(models, function(x) x$Ps)
    matrix_array = aperm(matrix_array, c(3, 4, 1, 2))
    dimnames(matrix_array)[[3]] = names(models)
  }
  Hs = llply(alply(matrix_array, 4, function(x) alply(x, 3)), function(x) KrzSubspace(x)$H)
  avgH = Reduce("+", Hs)/length(Hs)
  avgH.vec <- eigen(avgH)$vectors
  MCMC.H.val = laply(Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))

  # confidence intervals for variation in shared subspace directions
  observed = as.data.frame(HPDinterval(as.mcmc(MCMC.H.val), prob = prob))
  observed$mean = colMeans(MCMC.H.val)
  if(!is.null(ped)){
    n = dim(matrix_array)[1]
    m = dim(matrix_array)[3]
    MCMCsamp = dim(matrix_array)[4]
    traitnames = dimnames(matrix_array)[[1]]
    Gnames = dimnames(matrix_array)[[3]]
    rand.g_models <- array(NA ,c(n,n,m,MCMCsamp))
    dimnames(rand.g_models) <- list(traitnames, traitnames, Gnames)
    n_IDs = laply(IDs, length)
    filter_ped = function(x, i) x[rownames(x) %in% IDs[[i]],]
    for (i in 1:MCMCsamp){
      list_bv = vector("list", length(Gnames))
      names(list_bv) = Gnames
      for (line in Gnames){
        list_bv[[line]] <- filter_ped(rbv(ped, matrix_array[,,line,i]), line)
      }
      a.pop <- c(1, cumsum(n_IDs))
      pop.bv <- do.call(rbind, list_bv)
      rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
      for(j in seq_along(Gnames)){
        rand.g_models[,,j,i] <- cov(rand.pop.bv[a.pop[j]:a.pop[j+1],])
      }
    }
    null_Hs = llply(alply(rand.g_models, 4, function(x) alply(x, 3)), function(x) KrzSubspace(x)$H)
    null_avgH = Reduce("+", null_Hs)/length(null_Hs)
    null_avgH.vec <- eigen(null_avgH)$vectors
    null_MCMC.H.val = laply(null_Hs, function(mat) diag(t(avgH.vec) %*% mat %*% avgH.vec))
    null = as.data.frame(HPDinterval(as.mcmc(null_MCMC.H.val), prob = prob))
    null$mean = colMeans(null_MCMC.H.val)
    out = rbind(cbind(rank = 1:35, observed, type = "Observed"),
                cbind(rank = 1:35, null, type = "Randomized"))
  }
  else out = rbind(cbind(rank = 1:35, observed, type = "Observed"))
  return(out)
}
