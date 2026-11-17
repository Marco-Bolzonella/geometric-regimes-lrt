# Table S1, caption: largest discrepancy between the standard and
# the reference quadrature grids found by a dense search where the 30 x 25
# grid of quadrature_selfcheck.R finds its largest errors, at small xi and
# phi near 0 degrees (the law is symmetric under phi -> 180 - phi). The
# discrepancy oscillates on a fine scale, so the search ends with a local
# refinement around the ten best points.
#
# Usage: Rscript scripts/quadrature_finesearch.R [levels]   (default 1 2)

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))

args <- commandArgs(trailingOnly = TRUE)
LEVELS <- if (length(args)) as.integer(args) else 1:2
xi_grid <- 10^seq(log10(0.01), log10(3), length.out = 100)
phi_grid <- seq(1, 30, by = 0.25)

disc <- function(xi, phi_deg, w) {
  phi <- phi_deg * pi / 180
  100 * abs(tail_prob(w, xi, phi) - tail_prob(w, xi, phi, r_max = 12, n_theta = 192, n_r = 800))
}

rows <- list()
for (k in LEVELS) {
  w <- qchisq(2 * pnorm(k) - 1, 1)
  g <- expand.grid(xi = xi_grid, phi_deg = phi_grid)
  g$d <- mapply(disc, g$xi, g$phi_deg, MoreArgs = list(w = w))
  b <- g[which.max(g$d), ]
  # refine around the ten largest values of the dense grid
  starts <- g[order(-g$d), ][1:10, ]
  h <- do.call(rbind, lapply(seq_len(nrow(starts)), function(i)
    expand.grid(xi = starts$xi[i] * 10^seq(-0.05, 0.05, length.out = 41),
                phi_deg = starts$phi_deg[i] + seq(-1, 1, by = 0.05))))
  h <- h[h$phi_deg >= 1, ]
  h$d <- mapply(disc, h$xi, h$phi_deg, MoreArgs = list(w = w))
  bb <- h[which.max(h$d), ]
  cat(sprintf("%d sigma: dense grid max %.3f pp at xi = %.4f, phi = %.2f; refined max %.3f pp at xi = %.4f, phi = %.2f\n",
              k, b$d, b$xi, b$phi_deg, bb$d, bb$xi, bb$phi_deg))
  rows[[length(rows) + 1]] <- data.frame(sigma = k, dense_max_pp = b$d, dense_xi = b$xi,
                                         dense_phi_deg = b$phi_deg, refined_max_pp = bb$d,
                                         refined_xi = bb$xi, refined_phi_deg = bb$phi_deg)
}
out <- file.path(root, "tables", if (length(LEVELS) == 1) sprintf("quadrature_finesearch_%d.csv", LEVELS) else "quadrature_finesearch.csv")
write.csv(do.call(rbind, rows), out, row.names = FALSE)
cat(sprintf("wrote %s\n", out))
