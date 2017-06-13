#' Ratones
#'
#' Companion package for the Evolt... manuscript.
#'
#' @name ratones
#' @docType package
NULL

#' Linear distances for five mouse lines
#'
#' Skull distances measured from landmarks in 5 lines.
#'
#' @name ratones
#'
#' @docType data
#'
#' @usage data(ratones)
#'
#' @format List
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(ratones)
"ratones"

#' Linear distances for five mouse lines
#'
#' Skull distances measured from landmarks in 5 lines.
#'
#' @name ratonesdf
#'
#' @docType data
#'
#' @usage data(ratonesdf)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(ratonesdf)
#' delta_Z = colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "upwards",], IS_PM:BA_OPI)) -
#'           colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "downwards",], IS_PM:BA_OPI))
"ratonesdf"

#' Pedigree information for five mouse lines
#'
#' List with the pedigree for full experiment and relationship matrix for all measured animals
#'
#' @name ratones_ped
#'
#' @docType data
#'
#' @usage data(ratones_ped)
#'
#' @format list
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(ratones_ped)
#' if(!require(pedRSP)) devtools::install_github("jraffa/pedRSP");
#' library(pedRSP)
#' rsp = computeRSPs(kmat = ratones_ped$A, thre = 0.5)
"ratones_ped"

#' Posterior draws for BSFG models
#'
#' G and P matrices for the skull distances for all 5 lines obtained from the BSFG model.
#' The generateAnimalModelInput function will create the inputs for the MATLAB code available in
#' http://www2.stat.duke.edu/~sayan/bfgr/index.shtml it you wish to generate these yourself.
#'
#' @name ratones_models
#'
#' @docType data
#'
#' @usage data(ratones_models)
#'
#' @format list
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(ratones_models)
"ratones_models"

#' Posterior draws for BSFG within-groups model, including all 5 lines
#'
#' G and P matrices for the skull distances joining all 5 lines in single model for within-group matrices.
#' The generateAnimalModelInput function will create the inputs for the MATLAB code available in
#' http://www2.stat.duke.edu/~sayan/bfgr/index.shtml it you wish to generate these yourself.
#'
#' @name ratones_Wmodel
#'
#' @docType data
#'
#' @usage data(ratones_Wmodel)
#'
#' @format list
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(ratones_Wmodel)
"ratones_Wmodel"

#' Posterior draws for P matrix statistics
#'
#' P matrix statistics for all lines.
#'
#' @name P_stats
#'
#' @docType data
#'
#' @usage data(P_stats)
#'
#' @format list
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#'
#' @examples
#' data(P_stats)
#' \dontrun{
#' # To generate these from the posterior draws, run:
#' # (this can take a while)
#' # library(doMC)
#' # registerDoMC(5)
#' # P_stats = plyr::ldply(ratones_models, function(x) plyr::adply(x$Ps, 1,
#' # MeanMatrixStatistics, .progress = "text"), .parallel = TRUE)
#' # names(P_stats) <- gsub("pc1%", "pc1.percent", names(P_stats))
#' }
"P_stats"

#' Posterior draws for G matrix statistics
#'
#' G matrix statistics for all lines.
#'
#' @name G_stats
#'
#' @docType data
#'
#' @usage data(G_stats)
#'
#' @format list
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source \href{http://dryad_link.com}{Dryad Archive}
#' @examples
#' data(G_stats)
#' \dontrun{
#' # To generate these from the posterior draws, run:
#' # (this can take a while)
#' # library(doMC)
#' # registerDoMC(5)
#' # G_stats = plyr::ldply(ratones_models, function(x) plyr::adply(x$Gs, 1,
#' # MeanMatrixStatistics, .progress = "text"), .parallel = TRUE)
#' # names(G_stats) <- gsub("pc1%", "pc1.percent", names(G_stats))
#' }
"G_stats"
