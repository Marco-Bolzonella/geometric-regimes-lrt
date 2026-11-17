# Section 3.3, Section 5.3 and Appendix B.4: the closed-form law, Eq. (14),
# of W(c0) when sin(phi) = 0, checked against the signal-plane quadrature
# at phi = 179.9 and 179.99 degrees, and its 68.3% quantile at the values of xi(c0)
# found at the endpoints of the intervals in Section 5.3.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))

# distribution function of W(c0) on the half-line with endpoint at distance a = 1/(2 xi)
cdf_halfline <- function(w, xi) {
  a <- 1 / (2 * xi)
  ifelse(w <= a^2, 2 * pnorm(sqrt(w)) - 1, pnorm((w + a^2) / (2 * a)) - pnorm(-sqrt(w)))
}
quantile_halfline <- function(p, xi)
  uniroot(function(w) cdf_halfline(w, xi) - p, c(1e-12, 60), tol = 1e-10)$root

cover <- pnorm(1:3) - pnorm(-(1:3))
rows <- list()
for (phi_deg in c(179.9, 179.99)) {
  phi <- phi_deg * pi / 180
  for (xi in c(0.1, 0.3, 0.5, 1, 2, 5, 20, 100, 1000)) {
    for (k in 1:3) {
      w_hl <- quantile_halfline(cover[k], xi)
      # size, under the reference-grid quadrature, of the test based on the closed-form quantile
      size <- tail_prob(w_hl, xi, phi, r_max = 12, n_theta = 192, n_r = 800)
      rows[[length(rows) + 1]] <- data.frame(phi_deg = phi_deg, xi = xi, k_sigma = k,
                                             w_halfline = w_hl, nominal = 1 - cover[k],
                                             size = size, diff = size - (1 - cover[k]))
    }
  }
}
res <- do.call(rbind, rows)
print(res, digits = 4)
for (p in unique(res$phi_deg)) {
  r <- res[res$phi_deg == p, ]
  cat(sprintf("phi = %g deg: max |size - nominal| = %.4f over all xi, %.4f for xi >= 0.5\n",
              p, max(abs(r$diff)), max(abs(r$diff[r$xi >= 0.5]))))
}

# phi = 179.9 deg lies outside the range checked in Table S1, so check the
# 1-sigma sizes there by simulation as well
set.seed(20260909)
for (xi in c(0.1, 0.5)) {
  B <- 2e6; zu <- rnorm(B); zv <- rnorm(B)
  W <- compute_w_vectorized(zu, zv, xi, 179.9 * pi / 180)
  size <- mean(W > quantile_halfline(cover[1], xi))
  cat(sprintf("simulation, xi = %g, phi = 179.9 deg: 1-sigma size %.4f (se %.4f), nominal %.4f\n",
              xi, size, sqrt(size * (1 - size) / B), 1 - cover[1]))
}

# Section 5.3: xi(c0) at the endpoints of the 68.3% intervals of the six
# coefficients with xi(0) > 100 (tables/endpoint_accuracy.csv)
ep <- read.csv(file.path(root, "tables", "endpoint_accuracy.csv"))
xr <- range(ep$xi[ep$level_sigma == 1])
cat(sprintf("endpoints: xi between %.2f and %.2f; half-line 68.3%% quantile %.3f and %.3f\n",
            xr[1], xr[2], quantile_halfline(cover[1], xr[1]), quantile_halfline(cover[1], xr[2])))
write.csv(res, file.path(root, "tables", "halfline_check.csv"), row.names = FALSE)
