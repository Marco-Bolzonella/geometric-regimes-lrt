# Partial corrections of the chibar2 quantile (Appendix F): at both
# chibar2 endpoints of each 68.3% interval, the exact threshold against
# thresholds that adapt only xi(c0) or only phi(c0).

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R")); source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
raw <- load_raw(file.path(root, "data")); M <- whitening_matrix(raw$cov)
ALPHA1 <- 2 * pnorm(-1)
qW <- function(xi, phi) quantile_W(xi, phi, ALPHA1, r_max = 12, n_theta = 192, n_r = 800)
geo <- function(op, c0) {
  l <- sapply(raw$bins, function(b) raw$param[[b]][[op]]); q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
  T_ <- as.numeric(M %*% (l + 2 * c0 * q)); N_ <- as.numeric(2 * M %*% q); nT <- sqrt(sum(T_^2)); a <- sum(N_ * T_) / nT
  c(xi = sqrt(sum(N_^2)) / nT^2, phi = atan2(sqrt(sum((N_ - a * T_ / nT)^2)), a))
}
ends <- read.csv(file.path(root, "tables", "measured_intervals_1.csv"))
rows <- list()
for (op in ends$operator) {
  e <- ends[ends$operator == op, ]; g0 <- geo(op, 0)
  for (cb in c(e$lo_chibar, e$hi_chibar)) {
    g <- geo(op, cb)
    rows[[length(rows) + 1]] <- data.frame(operator = op, xi0 = g0["xi"], c0 = cb, xi = g["xi"], phi_deg = g["phi"] * 180 / pi,
      thr_exact = qW(g["xi"], g["phi"]), thr_xi_only = qW(g["xi"], g0["phi"]), thr_phi_only = qW(g0["xi"], g["phi"]),
      thr_chibar = qchisq(1 - 2 * ALPHA1, 1), row.names = NULL)
  }
  r <- rows[[length(rows)]]
}
res <- do.call(rbind, rows)
res$err_xi_only <- res$thr_xi_only - res$thr_exact; res$err_phi_only <- res$thr_phi_only - res$thr_exact
res$err_chibar <- res$thr_chibar - res$thr_exact
write.csv(res, file.path(root, "tables", "partial_correction_all.csv"), row.names = FALSE)
print(res[, c("operator", "xi0", "c0", "xi", "phi_deg", "thr_exact", "thr_xi_only", "thr_phi_only", "thr_chibar")], digits = 3, row.names = FALSE)
