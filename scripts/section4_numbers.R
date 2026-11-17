# Numbers in Sections 4 and 5 and Appendix F not produced by the other scripts.
# Prints to the console.
# (1) Appendix F, cost of Monte Carlo calibration: the margin by which W(c0) exceeds the
#     exact 1-sigma threshold between the two components of the cQu8
#     acceptance set, and the standard error of a threshold simulated with
#     8e4 draws there, {alpha(1-alpha)/n}^{1/2}/f(w), with f the exact
#     density (Eq. 14) at the reference-grid threshold.
# (2) Appendix F, Poisson counts: the sign of the 1-sigma size error of
#     chi2_1 against xi_P (needs tables/poisson_robustness.csv).
# (3) Appendix F, Poisson counts: agreement in sign with the predictions
#     of Figure 4 wherever these exceed 2.5 points.
# (4) Section 4.2: chi2_1 error rate at phi = 90 degrees, xi in [0.5, 1.5].

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "geom_tubular.R"))
source(file.path(root, "core", "load_operators.R"))

# (1) cQu8 margin
raw <- load_raw(file.path(root, "data"))
D <- solve(raw$cov)
M <- whitening_matrix(raw$cov)
hep <- yaml.load_file(file.path(root, "data", "hepdata_134011_table2_abs_ptthad.yaml"))
Y <- sapply(hep$dependent_variables[[1]]$values, function(v) v$value) /
  sapply(hep$dependent_variables[[2]]$values, function(v) v$value)
r0 <- Y - 1
op <- "cQu8"
l <- sapply(raw$bins, function(b) raw$param[[b]][[op]])
q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
geometry_at <- function(c0) {
  T_ <- as.numeric(M %*% (l + 2 * c0 * q))
  N_ <- as.numeric(2 * (M %*% q))
  normT <- sqrt(sum(T_^2))
  a <- sum(N_ * T_) / normT
  Nperp <- N_ - a * T_ / normT
  c(xi = sqrt(sum(N_^2)) / normT^2, phi = atan2(sqrt(sum(Nperp^2)), a))
}
dfun <- function(c) { r <- r0 - c * l - c^2 * q; sum(r * (D %*% r)) }
cs <- seq(-1, 1, by = 1e-4)
dmin <- optimize(dfun, cs[which.min(sapply(cs, dfun))] + c(-0.01, 0.01))$objective
alpha <- 1 - (pnorm(1) - pnorm(-1))
n_sim <- 8e4
c0s <- seq(-0.125, -0.065, by = 0.005)  # inside the excluded stretch [-0.13, -0.06]
res <- t(sapply(c0s, function(c0) {
  g <- geometry_at(c0)
  thr <- quantile_W(g["xi"], g["phi"], alpha, r_max = 12, n_theta = 192, n_r = 800)
  se <- sqrt(alpha * (1 - alpha) / n_sim) / fW_tubular(thr, g["xi"], g["phi"])
  c(c0 = c0, margin = dfun(c0) - dmin - thr, se = se)
}))
k <- which.max(res[, "margin"])
cat(sprintf("cQu8, 1 sigma: largest margin %.3f at c0 = %.3f, standard error of a threshold simulated with %.0e draws %.4f, ratio %.1f\n",
            res[k, "margin"], res[k, "c0"], n_sim, res[k, "se"], res[k, "margin"] / res[k, "se"]))
cat(sprintf("ratio at the grid points nearest the ends of the stretch: %.1f and %.1f\n",
            res[1, "margin"] / res[1, "se"], res[nrow(res), "margin"] / res[nrow(res), "se"]))

# (2) Poisson: sign of the chi2_1 size error at 1 sigma against xi_P
pr <- read.csv(file.path(root, "tables", "poisson_robustness.csv"))
p1 <- pr[pr$level_sigma == 1, ]
err <- p1$size_chi2_pct - p1$nominal_pct
cat(sprintf("Poisson, chi2_1 at 1 sigma: oversized for xi_P in [%.2f, %.2f], by up to %.1f points; undersized for xi_P in [%.2f, %.2f], by up to %.1f points\n",
            min(p1$xi[err > 0]), max(p1$xi[err > 0]), max(err),
            min(p1$xi[err < 0]), max(p1$xi[err < 0]), -min(err)))

# (3) Appendix F, Poisson counts: the sign of each fixed-threshold size
#     error agrees with the Gaussian-model prediction at (xi_P, phi_P)
#     wherever that prediction exceeds 2.5 points
for (k in 1:2) {
  s <- pr[pr$level_sigma == k, ]
  cn <- 2 * pnorm(k) - 1
  for (m in c("chi2", "chernoff")) {
    w <- if (m == "chi2") qchisq(cn, 1) else qchisq(2 * cn - 1, 1)
    pred <- vapply(seq_len(nrow(s)), function(i)
      100 * (tail_prob(w, s$xi[i], s$phi_deg[i] * pi / 180) - (1 - cn)), numeric(1))
    obs <- s[[paste0("size_", m, "_pct")]] - 100 * (1 - cn)
    big <- abs(pred) > 2.5
    cat(sprintf("Poisson, %d sigma, %s: predicted error above 2.5 points in %d rows, sign agrees in %d\n",
                k, m, sum(big), sum(sign(pred[big]) == sign(obs[big]))))
  }
}

# (4) Section 4.2: chi2_1 error rate at phi = 90 degrees for xi in [0.5, 1.5]
xis <- seq(0.5, 1.5, by = 0.05)
e90 <- vapply(xis, function(x) 100 * (1 - tail_prob(1, x, pi / 2, r_max = 12, n_theta = 192, n_r = 800) - (2 * pnorm(1) - 1)), numeric(1))
cat(sprintf("chi2_1 at 1 sigma, phi = 90 degrees, xi in [0.5, 1.5]: under-coverage %.2f to %.2f points\n",
            min(-e90), max(-e90)))
