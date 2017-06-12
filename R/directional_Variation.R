#' Directional variation
#'
#' Calculate the Evolvability and Conditional evolvability in the direction of selection.
#'
#' @param cov.matrix covariance matrix
#' @param line current line
#' @param delta_Z direction in phenotype space
#' @param Wmat optional fixed matrix for selection gradient reconstruction
#' @export
#' @importFrom evolqg ExtendMatrix Evolvability ConditionalEvolvability Normalize
#' @examples
#'
#' delta_Z = colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "upwards",], IS_PM:BA_OPI)) -
#'           colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "downwards",], IS_PM:BA_OPI))
#' \dontrun{
#' # this can take a while
#' library(doMC)
#' registerDoMC(5)
#' p_directional_stats <- ldply(ratones_models, function(model) adply(model$Ps, 1,
#'                                                                    directionalVariation,
#'                                                                    model$line,
#'                                                                    delta_Z), .parallel = TRUE)
#' DzPC1 = densityPlot(p_directional_stats, "DZpc1",
#'                     expression(paste("Vector correlation of ", delta, "z and E1")))
#' evolDZ = densityPlot(p_directional_stats, "evolDZ", "Scaled directional\n evolvability") +
#'            theme(legend.position = "none", text = element_text(size = 20))
#' condevolDZ = densityPlot(p_directional_stats, "condevolDZ",
#'                          "Scaled directional\n conditional evolvability") +
#'                theme(legend.position = "none", text = element_text(size = 20))
#' figure_4 <- ggdraw() + draw_plot(evolDZ, 0, 0.5, 0.5, 0.5) +
#'  draw_plot(condevolDZ, 0.5, 0.5, 0.5, 0.5) + draw_plot(DzPC1, 0.2, 0, 0.5, 0.5) +
#'  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0.2), c(1, 1, 0.5), size = 20)
#'  }
directionalVariation <- function(cov.matrix, line, delta_Z, Wmat = cov.matrix){
  beta_s <- solve(Wmat, delta_Z)
  beta_NN <- solve(ExtendMatrix(Wmat)[[1]], delta_Z)
  condEvol = (sum(ConditionalEvolvability(cov.matrix)))
  data.frame(evolDZ = Evolvability(cov.matrix, Normalize(delta_Z)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             evolBeta = Evolvability(cov.matrix, Normalize(beta_NN)) / (sum(diag(cov.matrix))/ncol(cov.matrix)),
             condevolDZ = ConditionalEvolvability(cov.matrix, Normalize(delta_Z)) / condEvol,
             condevolBeta = ConditionalEvolvability(cov.matrix, Normalize(beta_NN)) / condEvol,
             BetaPC1 = abs(vectorCor(beta_NN, eigen(cov.matrix)$vector[,1])),
             DZpc1 = abs(vectorCor(delta_Z, eigen(cov.matrix)$vector[,1])))
}
