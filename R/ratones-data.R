#' Ratones package
#'
#' Companion package for the "The Evolution of Phenotypic Integration: How directional selection reshapes covariation in mice" manuscript.
#'
#' @name ratones-package
#' @docType package
NULL

#' Data for five mouse lines
#'
#' Selection regime for body weight at 49 days was conducted in s and h lines for reduction, in s' and h' lines for increase, and a control line was maintained with random mating along generations, but avoiding full-sib matings. In downwards s and upwards s' lines inbreeding was performed by limiting population size and in downwards h and upwards h' full-sib mating was performed only during first generations.
#'
#' List object organized by line, containing individual information and data.
#'
#' \itemize{
#'  \item info.raw: complete information for all measured individuals, including both replicas.
#'          \itemize{
#'           \item ID: individual catalog number,
#'           \item SEX: sex (F: female, M: male),
#'           \item GER: generation of the long term experiment,
#'           \item LIN: experimental line (control t, downwards s, downwards h, upwards sp and upwards hp),
#'           \item MADRE: pedigree Dam,
#'           \item PADRE: pedigree Sire,
#'           \item AGE: age at sacrifice,
#'           \item P49: weight at 49 days of age,
#'           \item line: experimental line (control t, downwards s, downwards h, upwards s' and upwards h'),
#'           \item TAKE: paired individual replicas,
#'           \item selection: direction of selection (downwards, upwards or control),
#'           \item IS_PM-BA_OPI: set of 35 cranial euclidean distances, as indicated in the manuscript figure S4.
#'            }
#'  \item ed.raw: complete set of euclidean distances for all measured individuals, including replicas.
#'  \item info: information for all measured individuals.
#'  \item ed: individual mean between paired replicas of the 35 cranial euclidean distances.
#'  \item reps: vector of repeatabilities for the 35 cranial euclidean distances.
#'  \item model: linear model contoling for fixed effects (sex, age and generation).
#'  \item ed.means: vector of the line mean cranial euclidean distances.
#'  \item full: individual information, data and cranial euclinian distances'means between replicas.
#' \item gm_mean: line cranial geometric mean.
#' }
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
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#' @examples
#' data(ratones)
"ratones"

#' Information and cranial euclidean distances for five mouse lines
#'
#' Database for all individuals information and cranial euclidean distances mean between replicas.
#'
#' \itemize{
#'           \item .id: selection and line combined information.
#'           \item ID: individual catalog number,
#'           \item SEX: sex (F: female, M: male),
#'           \item GER: generation of the long term experiment,
#'           \item LIN: experimental line (control t, downwards s, downwards h, upwards sp and upwards hp),
#'           \item MADRE: pedigree Dam,
#'           \item PADRE: pedigree Sire,
#'           \item AGE: age at sacrifice,
#'           \item P49: weight at 49 days of age,
#'           \item line: experimental line,
#'           \item selection: direction of selection (downwards, upwards or control),
#'           \item IS_PM-BA_OPI: set of 35 cranial euclidean distances, as indicated in the manuscript figure S4.
#'            }
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
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#' @examples
#' data(ratonesdf)
#' delta_Z = colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "upwards",], IS_PM:BA_OPI)) -
#'           colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "downwards",], IS_PM:BA_OPI))
"ratonesdf"


#' Information and cranial euclidean distances for five mouse lines
#'
#' Database for all individuals information and cranial euclidean distances for both replicas.
#'
#' \itemize{
#'           \item .id: selection and line combined information.
#'           \item ID: individual catalog number,
#'           \item SEX: sex (F: female, M: male),
#'           \item GER: generation of the long term experiment,
#'           \item LIN: experimental line (control t, downwards s, downwards h, upwards sp and upwards hp),
#'           \item MADRE: pedigree Dam,
#'           \item PADRE: pedigree Sire,
#'           \item AGE: age at sacrifice,
#'           \item P49: weight at 49 days of age,
#'           \item line: experimental line,
#'           \item TAKE: paired individual replicas (1 and 2),
#'           \item selection: direction of selection (downwards, upwards or control),
#'           \item IS_PM-BA_OPI: set of 35 cranial euclidean distances, as indicated in the manuscript figure S4.
#'            }
#'
#' @name ratonesdf_raw
#'
#' @docType data
#'
#' @usage data(ratonesdf_raw)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#' @examples
#' data(ratonesdf)
#' delta_Z = colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "upwards",], IS_PM:BA_OPI)) -
#'           colMeans(dplyr::select(ratonesdf[ratonesdf$selection == "downwards",], IS_PM:BA_OPI))
"ratonesdf_raw"

#' Pedigree information for five mouse lines
#'
#' List with the pedigree for full experiment and relationship matrix for all measured animals.
#'
#' \itemize{
#'           \item ID: individual catalog number,
#'           \item MADRE: pedigree Dam,
#'           \item PADRE: pedigree Sire.
#'           }
#'
#' @name ratones_ped
#'
#' @docType data
#'
#' @usage data(ratones_ped)
#'
#' @format List
#'
#' @keywords datasets
#'
#' @references Penna et al. (2017) Evolution
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/**********}{PubMed})
#'
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#' @examples
#' data(ratones_ped)
#' if(!require(pedRSP)) devtools::install_github("jraffa/pedRSP");
#' library(pedRSP)
#' rsp = computeRSPs(kmat = ratones_ped$A, thre = 0.5)
"ratones_ped"

#' Posterior draws for BSFG models
#'
#' G and P matrices for the cranial distances for all 5 lines obtained from the BSFG model.
#' The generateAnimalModelInput function will create the inputs for the MATLAB code available in
#' http://www2.stat.duke.edu/~sayan/bfgr/index.shtml if you wish to generate these yourself.
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
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
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
#'
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#' @examples
#' data(ratones_Wmodel)
"ratones_Wmodel"

#' Posterior draws for P matrix evolutionary statistics
#'
#' P matrix evolutionary statistics for all lines.
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
#' Melo, D., G. Garcia, A. Hubbe, A. P. Assis, and G. Marroig. 2016. "EvolQG - An R Package for Evolutionary Quantitative Genetics [version 3; Referees: 2 Approved, 1 Approved with Reservations]." F1000Research 4: 925.
#'
#' Hansen, T. F., and Houle, D. (2008). Measuring and comparing evolvability and constraint in multivariate characters. Journal of evolutionary biology, 21(5), 1201-19. doi:10.1111/j.1420-9101.2008.01573.x
#'
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
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

#' Posterior draws for G matrix evolutionary statistics
#'
#' G matrix evolutionary statistics for all lines.
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
#' Melo, D., G. Garcia, A. Hubbe, A. P. Assis, and G. Marroig. 2016. "EvolQG - An R Package for Evolutionary Quantitative Genetics [version 3; Referees: 2 Approved, 1 Approved with Reservations]." F1000Research 4: 925.
#'
#' Hansen, T. F., and Houle, D. (2008). Measuring and comparing evolvability and constraint in multivariate characters. Journal of evolutionary biology, 21(5), 1201-19. doi:10.1111/j.1420-9101.2008.01573.x
#'
#' @source  Data from: The evolution of phenotypic integration: how directional selection reshapes covariation in mice. \href{http://dx.doi.org/10.5061/dryad.5gr8r}{Dryad Archive} doi:10.5061/dryad.5gr8r.
#'
#'
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
