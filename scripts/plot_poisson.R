# Figure 6, redrawn from tables/poisson_robustness.csv (written by
# make_poisson_robustness.R).

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "plot_theme.R"))

MU_SCALES <- c(1, 3, 10, 30, 100)
result_df <- read.csv(file.path(root, "tables", "poisson_robustness.csv"))
plot_df <- do.call(rbind, lapply(c(exact = "size_exact_pct", chi2 = "size_chi2_pct",
                                   chibar = "size_chernoff_pct"), function(col) {
  data.frame(result_df[, c("operator", "regime", "level_sigma", "nominal_pct", "mu_scale")],
             rate = result_df[[col]], method = col)
}))
plot_df$method <- factor(plot_df$method, levels = c("size_exact_pct", "size_chi2_pct", "size_chernoff_pct"),
                         labels = c("'Exact quantile'", "chi[1]^2~quantile", "bar(chi)^2~quantile"))
plot_df$level <- factor(plot_df$level_sigma, levels = 1:2,
                        labels = c("'1'*sigma~'(nominal 31.7%)'", "'2'*sigma~'(nominal 4.6%)'"))
plot_df$regime <- factor(paste("Regime", plot_df$regime), levels = paste("Regime", c("I", "II", "III")))
nom_df <- unique(plot_df[, c("level", "nominal_pct")])
p <- ggplot(plot_df, aes(mu_scale, rate, group = operator, colour = regime)) +
  geom_hline(data = nom_df, aes(yintercept = nominal_pct), linetype = "dashed", colour = "grey30") +
  geom_line(linewidth = 0.35) + geom_point(size = 0.7) +
  scale_x_log10(breaks = MU_SCALES) +
  scale_colour_manual(values = COL_REGIME, name = NULL) +
  facet_grid(level ~ method, scales = "free_y", labeller = label_parsed) +
  labs(x = "Expected count per bin under the null, m (log scale)", y = "Error rate (%)") +
  theme_paper() + theme(legend.position = "bottom")
save_fig(p, file.path(root, "figures", "poisson_size.png"), height = 4.2)
