# Figure 7, Table 1 (MLE and W(0)) and Section 5.3: confidence sets for the 15 coefficients from the
# measured data (HEPData record 10.17182/hepdata.134011.v2, Table 2).
# Y_i = measured / Powheg+Pythia8 prediction, so mu_i^SM = 1. W(c0) comes
# from the 8-bin likelihood; the exact threshold is recomputed from
# (xi(c0), phi(c0)) at every c0. Crossings are found on a grid and refined
# with uniroot, so sets with several components are found too.
#
# Usage: Rscript scripts/measured_intervals.R [levels]   (default 1 2)

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))

args <- commandArgs(trailingOnly = TRUE)
LEVELS <- if (length(args)) as.integer(args) else c(1L, 2L)
N_GRID <- 300

load_ratio_data <- function(data_dir_path) {
  d <- yaml.load_file(file.path(data_dir_path, "hepdata_134011_table2_abs_ptthad.yaml"))
  meas <- sapply(d$dependent_variables[[1]]$values, function(v) v$value)
  pred <- sapply(d$dependent_variables[[2]]$values, function(v) v$value)
  stopifnot(grepl("PWG\\+PY8", d$dependent_variables[[2]]$header$name))
  meas / pred
}

raw <- load_raw(file.path(root, "data"))
D <- solve(raw$cov)
M <- whitening_matrix(raw$cov)
Y <- load_ratio_data(file.path(root, "data"))
r0 <- Y - 1
chi2_sm <- as.numeric(t(r0) %*% D %*% r0)
cat(sprintf("data ratio: %s\nSM chi2 = %.2f on %d df\n",
            paste(sprintf("%.6f", Y), collapse = " "), chi2_sm, length(Y)))

geometry_at <- function(c0, l, q) {
  T_ <- as.numeric(M %*% (l + 2 * c0 * q))
  N_ <- as.numeric(2 * (M %*% q))
  normT <- sqrt(sum(T_^2))
  a <- sum(N_ * T_) / normT
  Nperp <- N_ - a * T_ / normT
  c(xi = sqrt(sum(N_^2)) / normT^2, phi = atan2(sqrt(sum(Nperp^2)), a))
}

make_W <- function(l, q) {
  # d(c) = (r0 - c l - c^2 q)' D (r0 - c l - c^2 q)
  dl <- as.numeric(D %*% l); dq <- as.numeric(D %*% q)
  k0 <- sum(r0 * D %*% r0); k1 <- -2 * sum(r0 * dl)
  k2 <- sum(l * dl) - 2 * sum(r0 * dq); k3 <- 2 * sum(l * dq); k4 <- sum(q * dq)
  dfun <- function(c) k0 + k1 * c + k2 * c^2 + k3 * c^3 + k4 * c^4
  st <- polyroot(c(k1, 2 * k2, 3 * k3, 4 * k4))
  st <- Re(st[abs(Im(st)) < 1e-8])
  chat <- st[which.min(dfun(st))]
  list(W = function(c0) dfun(c0) - dfun(chat), chat = chat)
}

# reference grid
REF <- list(r_max = 12, n_theta = 192, n_r = 800)
thr_exact <- function(c0, l, q, alpha) {
  g <- geometry_at(c0, l, q)
  quantile_W(g["xi"], g["phi"], alpha, r_max = REF$r_max, n_theta = REF$n_theta, n_r = REF$n_r)
}

# acceptance set as a matrix of [lo, hi] rows
acceptance <- function(Wf, thr, lo, hi, n = N_GRID) {
  grid <- seq(lo, hi, length.out = n)
  g <- vapply(grid, function(c) Wf(c) - thr(c), numeric(1))
  inside <- g <= 0
  stopifnot(!inside[1], !inside[n])
  cross <- which(diff(inside) != 0)
  pts <- vapply(cross, function(k) {
    uniroot(function(c) Wf(c) - thr(c), grid[k + 0:1], tol = 1e-7)$root
  }, numeric(1))
  matrix(pts, ncol = 2, byrow = TRUE)
}

ops <- build_table1(file.path(root, "data"))
rows <- list()
for (lev in LEVELS) {
  cover <- pnorm(lev) - pnorm(-lev)
  alpha <- 1 - cover
  w_chi2 <- qchisq(cover, 1)
  w_bar <- qchisq(max(2 * cover - 1, 0), 1)
  for (op in ops$operator) {
    l <- sapply(raw$bins, function(b) raw$param[[b]][[op]])
    q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
    wl <- make_W(l, q)
    # scan until W is well above every threshold on both sides
    big <- 3 * w_chi2 + 2
    step <- 0.05
    lo <- wl$chat - step; while (wl$W(lo) < big) { lo <- lo - step; step <- step * 1.5 }
    step <- 0.05
    hi <- wl$chat + step; while (wl$W(hi) < big) { hi <- hi + step; step <- step * 1.5 }
    ex <- acceptance(wl$W, function(c) thr_exact(c, l, q, alpha), lo, hi)
    ch <- acceptance(wl$W, function(c) w_chi2, lo, hi, n = 4000)
    cb <- acceptance(wl$W, function(c) w_bar, lo, hi, n = 4000)
    main <- function(s) { k <- which.max(s[, 2] - s[, 1]); s[k, ] }
    fmt_set <- function(s) paste(sprintf("[%.4f,%.4f]", s[, 1], s[, 2]), collapse = "U")
    total <- function(s) sum(s[, 2] - s[, 1])
    w_ex <- total(ex); w_ch <- total(ch); w_cb <- total(cb)
    g0 <- geometry_at(0, l, q)
    rows[[length(rows) + 1]] <- data.frame(
      level_sigma = lev, operator = op, xi0 = g0["xi"], phi0_deg = g0["phi"] * 180 / pi,
      chat = wl$chat, W0 = wl$W(0),
      width_exact = w_ex, width_chi2 = w_ch, width_chibar = w_cb,
      pct_chi2 = 100 * (w_ch / w_ex - 1), pct_chibar = 100 * (w_cb / w_ex - 1),
      lo_exact = main(ex)[1], hi_exact = main(ex)[2],
      lo_chibar = main(cb)[1], hi_chibar = main(cb)[2],
      n_components = nrow(ex), set_exact = fmt_set(ex), set_chi2 = fmt_set(ch),
      set_chibar = fmt_set(cb), row.names = NULL)
    cat(sprintf("%dsigma %-6s xi0=%9.2f chat=%+.3f W0=%.3f exact=%.4f chi2=%+.1f%% chibar=%+.1f%%  %s\n",
                lev, op, g0["xi"], wl$chat, wl$W(0), w_ex, 100 * (w_ch / w_ex - 1),
                100 * (w_cb / w_ex - 1), fmt_set(ex)))
  }
}
res <- do.call(rbind, rows)
for (lev in LEVELS) {
  out <- file.path(root, "tables", sprintf("measured_intervals_%d.csv", lev))
  write.csv(res[res$level_sigma == lev, ], out, row.names = FALSE)
  cat(sprintf("wrote %s\n", out))
}
