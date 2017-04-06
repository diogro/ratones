source('./R/read_ratones.R')

full_data$gm
m_data = melt(select(full_data, P49, ID, LIN, SEX, IS_PM:BA_OPI, gm), id.vars = c("ID", "P49", "LIN", "SEX", "gm"))
current_var = "IS_PM"
llply(unique(m_data$variable), function(current_var){
all_regression = ggplot(filter(m_data, variable == current_var), aes(gm, value)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~LIN+SEX, scale = "free_y", ncol = 5)
save_plot(paste0("regression_", current_var, ".png"), all_regression, 
          base_height = 5, base_aspect_ratio = 2)})
