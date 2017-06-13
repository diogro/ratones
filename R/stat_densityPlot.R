#' Generate denisty plot for some matrix statistic
#'
#'@param mcmc_stats list of outputs from MeanMatrixStatistics
#'@param variable variable to be ploted, any column in the mcmc_stats objects
#'@param xlab x axis label
#'@return ggplot2 density plot object
#'@export
#'@importFrom ggplot2 ggplot geom_density aes aes_string scale_fill_manual labs theme element_text
#'@importFrom cowplot background_grid ggdraw draw_plot panel_border
#'@importFrom viridis viridis
#'@importFrom dplyr select_
#'@importFrom tidyr separate
#' @examples
#' data(P_stats)
#' library(ggplot2)
#' library(cowplot)
#' r2_plot = densityPlot(P_stats, "MeanSquaredCorrelation",
#'                       "Mean squared correlation")
#' pc1.percent_plot = densityPlot(P_stats, "pc1.percent",
#'                                "Proportion of variation in E1") +
#'                      theme(legend.position = "none")
#' flexibility_plot = densityPlot(P_stats, "flexibility",
#'                                "Mean flexibility") +
#'                      theme(legend.position = "none")
#' evolvability_plot = densityPlot(P_stats, "evolvability",
#'                                 "Mean evolvability") +
#'                       theme(legend.position = "none")
#' cond_evolvability_plot = densityPlot(P_stats, "conditional.evolvability",
#'                                      "Mean conditional evolvability") +
#'                            theme(legend.position = "none")
#' figure_3 <- ggdraw() +
#'     draw_plot(r2_plot, 0, 0.5, 0.5, 0.5) +
#'     draw_plot(pc1.percent_plot, 0.5, 0.5, 0.5, 0.5) +
#'     draw_plot(flexibility_plot, 0, 0, 0.5, 0.5) +
#'     draw_plot(evolvability_plot, 0.5, 0, 0.5, 0.5) +
#'     draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5),
#'                     c(1, 1, 0.5, 0.5), size = 20)
densityPlot <- function(mcmc_stats, variable, xlab){
  global_stats <- separate(dplyr::select_(mcmc_stats, ".id", variable), ".id", c('selection', 'line'), sep = "\\.")
  lines = c("t", "h'", "s'", "h", "s")
  global_stats$line <- factor(global_stats$line, levels = lines)
  selection = line = NULL
  plot <- ggplot(global_stats, aes_string(variable)) +
    geom_density(alpha = 0.5, aes(group = interaction(selection, line),
                                  fill = interaction(selection, line))) +
    scale_fill_manual(values = viridis(5), name = "Line",
                      labels = c("Control t", "Upwards h'", "Upwards s'", "Downwards h", "Downwards s")) +
    background_grid(major = 'x', minor = "none") +
    panel_border() +
    labs(x = xlab, y = "Density") +
    theme(legend.position = c(0.7, 0.7), text = element_text(size = 20))
}
