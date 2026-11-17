# Figure 5 (main text: exact quantile only) and Figure S1
# (Supplementary Material: exact, chi2_1 and chibar2 quantiles), redrawn
# from tables/robustness_grid_levels.rds, written by make_robustness_heatmaps.R.
# Distributions in columns, quantiles in rows.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "plot_theme.R"))

g <- readRDS(file.path(root, "tables", "robustness_grid_levels.rds"))
operators <- build_table1(file.path(root, "data"))
scenarios <- names(g$results[[1]])
scenario_labels <- c(Normal = "Normal", t5 = "t[5]", t3 = "t[3]")
method_labels <- c(exact = "'Exact quantile'", chi2 = "chi[1]^2~quantile", chibar = "bar(chi)^2~quantile")

k2 <- which(g$levels == 2)
nom <- 100 * g$alphas[k2]
dev <- 4
df <- do.call(rbind, lapply(names(method_labels), function(m) do.call(rbind, lapply(scenarios, function(s) {
  d <- expand.grid(xi = g$xi_grid, phi_deg = g$phi_grid)
  d$rate <- as.vector(t(g$results[[k2]][[s]][[m]])) * 100
  d$scenario <- scenario_labels[[s]]
  d$method <- method_labels[[m]]
  d
}))))
df$scenario <- factor(df$scenario, levels = unname(scenario_labels))
df$method <- factor(df$method, levels = unname(method_labels))

heatmap <- function(d)
  ggplot(d, aes(xi, phi_deg, fill = rate)) +
    geom_raster() +
    scale_x_xi() + scale_y_phi() +
    scale_fill_gradient2(name = "Error\nrate (%)", low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = nom, limits = c(max(0, nom - dev), nom + dev),
                         oob = scales::squish) +
    geom_operators(operators) +
    labs(x = LAB_XI, y = LAB_PHI) +
    theme_heatmap()

p_all <- heatmap(df) + facet_grid(method ~ scenario, labeller = label_parsed)
save_fig(p_all, file.path(root, "figures", "robustness_heatmaps_all.png"), height = 6.0)
p_exact <- heatmap(df[df$method == method_labels[["exact"]], ]) + facet_grid(~scenario, labeller = label_parsed)
save_fig(p_exact, file.path(root, "figures", "robustness_heatmaps.png"), height = 2.5)
