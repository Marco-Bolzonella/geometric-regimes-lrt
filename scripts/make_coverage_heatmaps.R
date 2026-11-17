# Figure 4 (Section 4.2): ratio of actual to nominal error rate of the tests based on the chi2_1 and chibar2 quantiles at
# 1 and 2 sigma over the (xi, phi) plane; 3-sigma summaries are printed.
# The grids are cached in tables/coverage_grids_levels.rds.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "plot_theme.R"))

N_XI  <- 80
N_PHI <- 70
LEVELS <- 1:3                                   # k sigma; 1 and 2 are plotted, 3 is summarised
xi_grid  <- 10^seq(log10(0.05), log10(30000), length.out = N_XI)
phi_grid <- seq(2, 178, length.out = N_PHI)

# coverage of a fixed threshold over the grid
achieved_grid <- function(w_threshold) {
  g <- expand.grid(xi = xi_grid, phi_deg = phi_grid)
  g$achieved <- vapply(seq_len(nrow(g)), function(i) {
    1 - tail_prob(w_threshold, g$xi[i], g$phi_deg[i] * pi / 180)
  }, numeric(1))
  g
}

cache <- file.path(root, "tables", "coverage_grids_levels.rds")
if (file.exists(cache)) {
  grids <- readRDS(cache)
} else {
  grids <- list()
  for (k in LEVELS) {
    cov_nom <- 2 * pnorm(k) - 1
    cat(sprintf("computing %d-sigma grids...\n", k))
    grids[[k]] <- list(cov_nom = cov_nom,
                       chi2 = achieved_grid(qchisq(cov_nom, df = 1)),
                       chernoff = achieved_grid(qchisq(2 * cov_nom - 1, df = 1)))
  }
  saveRDS(grids, cache)
}
for (k in LEVELS) {
  cn <- grids[[k]]$cov_nom
  for (m in c("chi2", "chernoff")) {
    a <- grids[[k]][[m]]$achieved
    pp <- 100 * (a - cn)
    rel <- 100 * ((1 - a) / (1 - cn) - 1)       # relative error of the tail probability (size)
    cat(sprintf("%d sigma %-8s: coverage error [%.3f, %.3f] pp; size error [%.1f, %.1f] %%\n",
                k, m, min(pp), max(pp), min(rel), max(rel)))
  }
}

operators <- build_table1(file.path(root, "data"))

# ratio of the actual to the nominal error rate of the test that uses a
# fixed quantile, on a log scale so that 1/2 and 2 are equally far from 1;
# one colour scale for all panels, laid out as Figure 6 (levels in rows,
# quantiles in columns)
df <- do.call(rbind, lapply(1:2, function(k) do.call(rbind, lapply(c("chi2", "chernoff"), function(m) {
  d <- grids[[k]][[m]]
  d$ratio <- (1 - d$achieved) / (1 - grids[[k]]$cov_nom)
  d$level <- sprintf("'%d'*sigma~'(nominal %.1f%%)'", k, 100 * (1 - grids[[k]]$cov_nom))
  d$method <- m
  d
}))))
df$method <- factor(df$method, levels = c("chi2", "chernoff"),
                    labels = c("chi[1]^2~quantile", "bar(chi)^2~quantile"))
p <- ggplot(df, aes(xi, phi_deg, fill = ratio)) +
  geom_raster() +
  scale_x_xi() + scale_y_phi() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       trans = "log2", limits = c(1 / 4, 4), oob = scales::squish,
                       breaks = c(1 / 4, 1 / 2, 1, 2, 4), labels = c("1/4", "1/2", "1", "2", "4"),
                       name = "Actual /\nnominal\nerror rate") +
  geom_operators(operators) +
  facet_grid(level ~ method, labeller = label_parsed) +
  labs(x = LAB_XI, y = LAB_PHI) +
  theme_heatmap()

save_fig(p, file.path(root, "figures", "error_rate_heatmaps.png"), height = 4.4)
