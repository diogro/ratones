load("./R/read_ratones.R")

m_full_data = melt(full_data, id.vars = names(full_data)[c(1:8, 10:12)])

full_trait_plots = ggplot(m_full_data %>% filter(variable != 'P49'), aes(original_line, value, group = original_line, fill = selection)) + geom_violin() + scale_fill_manual(values = c(c, dw, up)) + facet_wrap(~variable, scale = "free_y", ncol = 5) + labs(y = "Linear distance between landmarks (mm)", x = "line") + geom_jitter(height = 0, width = 0.08, alpha = 0.3)

#figure_folder = "~/Dropbox/labbio/Shared Lab/Ratones_shared/"
figure_folder = "~/Desktop/"

save_plot(paste0(figure_folder, "figureS4.pdf"), full_trait_plots, ncol = 5, nrow = 7, base_height = 3)

p49_full_data = m_full_data %>% filter(variable == 'P49')
p49_full_data$SEX %<>% {gsub("M", "Male", .)} %>% {gsub("F", "Female", .)}
p49_plot = ggplot(p49_full_data, aes(line, value, group = original_line, fill = selection)) + geom_violin() + scale_color_manual(values = c(c, dw, up)) + scale_fill_manual(values = c(c, dw, up)) + facet_wrap(~SEX) + background_grid(major = 'y', minor = "none") + labs(y = "Weight at 49 days (g)", x = "line") + geom_jitter(height = 0, width = 0.08, alpha = 0.3)

save_plot(paste0(figure_folder, "figureS2.png"), p49_plot, base_height = 4, base_aspect_ratio = 1.7)


traits = full_data %>% dplyr::select(IS_PM:BA_OPI)
full_data$gm = apply(traits, 1, gm_mean)
# p49_gm_plot = ggplot(full_data, aes(log(P49), gm, group = .id, color = .id)) + geom_point() + geom_smooth(method = "lm")
# save_plot("~/Dropbox/labbio/Shared Lab/Ratones_shared/p49_gm.pdf", p49_gm_plot, base_aspect_ratio = 1.6, base_height = 12)

gm_full_data = melt(full_data, id.vars = names(full_data)[c(1:8, 10:12)]) %>% filter(variable == 'gm')
gm_full_data$SEX %<>% {gsub("M", "Male", .)} %>% {gsub("F", "Female", .)}
gm_plot = ggplot(gm_full_data, aes(original_line, value, group = original_line, fill = selection)) + geom_violin() + scale_fill_manual(values = c(c, dw, up)) + facet_wrap(~SEX) + background_grid(major = 'y', minor = "none") + labs(y = "Geometric mean of cranial traits", x = "line")+ geom_jitter(height = 0, width = 0.08, alpha = 0.3)

save_plot(paste0(figure_folder, "figureS5.pdf"), gm_plot, base_height = 4, base_aspect_ratio = 1.7)

full_data %>% count(line, selection, SEX) %>% xtable
full_data %>% count(line, selection, GER) %>% xtable
full_data %>% count(line, selection) %>% xtable

cvs = data.frame(t(laply(names(main.data), function(x) sqrt(diag(r_models[[x]]$MAP))/ main.data[[x]]$ed.means)))
names(cvs) <- names(main.data)
cvs$traits = factor(rownames(cvs), levels = rownames(cvs)[order(cvs$control.t)])

cv_plot = ggplot(melt(cvs), aes(traits, value, group = variable, color = variable, shape = variable)) + geom_point() + geom_line() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Traits", y = "Coeficient of variation") + scale_color_manual(name = "line", values = viridis(5)) + scale_y_continuous(limits = c(0, 0.2)) + background_grid(major = 'y', minor = "y") + theme(legend.position = c(0.15, 0.8))

save_plot(paste0(figure_folder, "figureS6.pdf"), cv_plot, base_height = 5.5, base_aspect_ratio = 2)
