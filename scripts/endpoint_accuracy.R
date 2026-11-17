# Accuracy of the reference grid at the endpoints of the confidence
# intervals, where phi is close to 0 or 180 degrees: the coverage of the
# reference-grid threshold is recomputed from Eq. (15).

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R")); source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "geom_tubular.R")); source(file.path(root, "core", "load_operators.R"))
raw <- load_raw(file.path(root, "data")); M <- whitening_matrix(raw$cov)
geo <- function(op, c0) {
  l <- sapply(raw$bins, function(b) raw$param[[b]][[op]]); q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
  T_ <- as.numeric(M %*% (l + 2 * c0 * q)); N_ <- as.numeric(2 * M %*% q); nT <- sqrt(sum(T_^2)); a <- sum(N_ * T_) / nT
  c(xi = sqrt(sum(N_^2)) / nT^2, phi = atan2(sqrt(sum((N_ - a * T_ / nT)^2)), a))
}
rows <- list()
for (k in 1:2) {
  d <- read.csv(file.path(root, "tables", sprintf("measured_intervals_%d.csv", k))); cv <- 2 * pnorm(k) - 1
  for (op in c("cQj11", "ctj1", "cQu1", "cQd1", "cG", "ctu1")) {
    r <- d[d$operator == op, ]
    for (c0 in c(r$lo_exact, r$hi_exact)) {
      g <- geo(op, c0)
      w <- quantile_W(g["xi"], g["phi"], 1 - cv, r_max = 12, n_theta = 192, n_r = 800)
      tt <- tail_tubular(w, g["xi"], g["phi"])
      rows[[length(rows) + 1]] <- data.frame(level_sigma = k, operator = op, c0 = c0, xi = g["xi"],
        phi_deg = g["phi"] * 180 / pi, w_ref = w, cover_tubular = 1 - tt, diff_pp = 100 * ((1 - tt) - cv), row.names = NULL)
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(root, "tables", "endpoint_accuracy.csv"), row.names = FALSE)
cat(sprintf("phi within %.2f deg of 0 or 180; max |reference grid - Eq. (15)| = %.4f pp\n",
            max(pmin(res$phi_deg, 180 - res$phi_deg)), max(abs(res$diff_pp))))
