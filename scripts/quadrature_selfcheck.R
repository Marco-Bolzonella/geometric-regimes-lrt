# Table S1 (Section 4.1, Appendix F): standard grid against reference grid at 1-4 sigma.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))

N_XI <- 30
N_PHI <- 25
xi_grid <- 10^seq(log10(0.01), log10(1e5), length.out = N_XI)
phi_grid <- seq(1, 179, length.out = N_PHI) * pi / 180

sigma_levels <- 1:4
cov_nom <- 2 * pnorm(sigma_levels) - 1
w_thr <- qchisq(cov_nom, df = 1)

max_err <- numeric(length(sigma_levels))
worst_at <- vector("list", length(sigma_levels))

for (k in seq_along(sigma_levels)) {
  worst <- 0; worst_xi <- NA; worst_phi <- NA
  for (xi in xi_grid) {
    for (phi in phi_grid) {
      std  <- 1 - tail_prob(w_thr[k], xi, phi)
      conv <- 1 - tail_prob(w_thr[k], xi, phi, n_theta = 192, n_r = 800, r_max = 12)
      err <- 100 * abs(std - conv)
      if (err > worst) { worst <- err; worst_xi <- xi; worst_phi <- phi * 180 / pi }
    }
  }
  max_err[k] <- worst
  worst_at[[k]] <- c(xi = worst_xi, phi_deg = worst_phi)
  cat(sprintf("sigma=%d: max err = %.4f pp at xi=%.3g, phi=%.1f deg\n", sigma_levels[k], worst, worst_xi, worst_phi))
}

out <- data.frame(sigma = sigma_levels, cov_nom = cov_nom, max_abs_error_pp = max_err,
                   worst_xi = sapply(worst_at, `[`, "xi"), worst_phi_deg = sapply(worst_at, `[`, "phi_deg"))
write.csv(out, file.path(root, "tables", "quadrature_selfcheck.csv"), row.names = FALSE)
cat("DONE\n")
