# Table S2 and Figure 6 (Section 4.4): error rates of the tests based on the three quantiles when the
# data are Poisson counts, for the 15 coefficients and m = 1, 3, 10, 30, 100
# expected counts per bin. Thresholds come from the Gaussian model fitted to
# the counts (covariance I/m for Y/m); W(0) from the Gaussian working
# likelihood, the statistic the method calibrates.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "geom_poisson.R"))
library(ggplot2)

SEED <- 20260909
B <- 50000
LEVELS <- c(1, 2)
ALPHAS <- 2 * pnorm(-LEVELS)
MU_SCALES <- c(1, 3, 10, 30, 100)
C_GRID <- seq(-30, 30, length.out = 3000)

raw <- load_raw(file.path(root, "data"))
ops <- build_table1(file.path(root, "data"))
rownames(ops) <- ops$operator

# equivalent counts per bin of the measured data, m = 1/delta^2, with delta
# the relative statistical uncertainty (Section 4.4, last sentence)
hep <- yaml.load_file(file.path(root, "data", "hepdata_134011_table2_abs_ptthad.yaml"))
delta <- sapply(hep$dependent_variables[[1]]$values, function(v) {
  e <- Filter(function(x) identical(x$label, "stat"), v$errors)[[1]]
  as.numeric(sub("%", "", e$symerror)) / 100
})
cat("Equivalent counts per bin, 1/delta^2:", sprintf("%.0f", 1 / delta^2), "\n")

rows <- list()
seed_ctr <- SEED
for (op in ops$operator) {
  l <- sapply(raw$bins, function(b) raw$param[[b]][[op]])
  q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
  for (mu in MU_SCALES) {
    # geometry under Sigma_P = I/m
    T0 <- sqrt(mu) * l
    N  <- 2 * sqrt(mu) * q
    e1 <- T0 / sqrt(sum(T0^2))
    a  <- sum(N * e1)
    xi  <- sqrt(sum(N^2)) / sum(T0^2)
    phi <- atan2(sqrt(sum((N - a * e1)^2)), a)
    seed_ctr <- seed_ctr + 1
    res <- simulate_poisson_W(l, q, mu_scale = mu, B = B, c_grid = C_GRID, seed = seed_ctr)
    if (res$boundary_frac > 0.001) {
      cat(sprintf("WARNING: %s mu=%d boundary_frac=%.4f -- widen C_GRID\n", op, mu, res$boundary_frac))
    }
    for (k in seq_along(LEVELS)) {
      alpha <- ALPHAS[k]
      w_exact <- quantile_W(xi, phi, alpha)
      w_chi2 <- qchisq(1 - alpha, df = 1)
      w_chibar <- qchisq(1 - 2 * alpha, df = 1)
      rows[[length(rows) + 1]] <- data.frame(
        operator = op, regime = ops[op, "regime"], level_sigma = LEVELS[k], nominal_pct = 100 * alpha,
        mu_scale = mu, xi = xi, phi_deg = phi * 180 / pi,
        w_exact_threshold = w_exact,
        size_exact_pct = 100 * mean(res$W > w_exact),
        size_chi2_pct = 100 * mean(res$W > w_chi2),
        size_chernoff_pct = 100 * mean(res$W > w_chibar),
        se_pct = 100 * sqrt(alpha * (1 - alpha) / B), boundary_frac = res$boundary_frac)
    }
    r <- rows[[length(rows)]]
    cat(sprintf("%-6s mu=%3d xi_P=%10.2f  2sigma: exact=%.2f%% chi2_1=%.2f%% chibar2=%.2f%%\n",
                op, mu, xi, r$size_exact_pct, r$size_chi2_pct, r$size_chernoff_pct))
  }
}
result_df <- do.call(rbind, rows)
write.csv(result_df, file.path(root, "tables", "poisson_robustness.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\n", file.path(root, "tables", "poisson_robustness.csv")))
for (k in LEVELS) for (m in c("size_exact_pct", "size_chi2_pct", "size_chernoff_pct")) {
  v <- result_df[result_df$level_sigma == k, m]
  cat(sprintf("%d sigma %-18s range [%.2f, %.2f] %%\n", k, m, min(v), max(v)))
}

source(file.path(root, "scripts", "plot_poisson.R"))   # Figure 6
cat(sprintf("seed base=%d, B=%d per (coefficient, m)\n", SEED, B))
