# Coverage and length of the confidence sets (Section 5.4, Appendix F). For
# each coefficient and true value c in {0, chat}, B datasets are drawn from
# N(mu(c), Sigma0); for each one the acceptance set {c0 : W(c0) <= w_k(c0)}
# is built on a grid of c0 and checked for c. The exact quantile w_k depends
# on c0 only, so it is computed once per coefficient (reference grid) and
# interpolated. Coverage is also computed directly as P{W(c) <= w_k(c)}. For
# each method the script also records the mean and median length of the
# acceptance set, the mean offset of its centre from the true value, and the
# fraction of sets with more than one component.
#
# Usage: Rscript scripts/coverage_length.R [B] [quick]   (default B = 50000;
# "quick" runs two coefficients on coarse grids, as a smoke test)

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 50000L
QUICK <- length(args) > 1 && args[2] == "quick"   # smoke test: 2 coefficients, coarse grids
SEED <- 20260909
N_C <- 400
LEVELS <- c(1, 2)
REF <- list(r_max = 12, n_theta = 192, n_r = 800)

raw <- load_raw(file.path(root, "data"))
S0 <- raw$cov; D <- solve(S0); M <- whitening_matrix(raw$cov); Lc <- t(chol(S0))
ints <- lapply(1:2, function(k) read.csv(file.path(root, "tables", sprintf("measured_intervals_%d.csv", k))))

geometry_at <- function(c0, l, q) {
  T_ <- as.numeric(M %*% (l + 2 * c0 * q)); N_ <- as.numeric(2 * (M %*% q))
  nT <- sqrt(sum(T_^2)); a <- sum(N_ * T_) / nT
  c(xi = sqrt(sum(N_^2)) / nT^2, phi = atan2(sqrt(sum((N_ - a * T_ / nT)^2)), a))
}

# real roots of a3 x^3 + a2 x^2 + a1 x + a0 (a1, a0 vectors); NA if missing
cubic_roots_vec <- function(a3, a2, a1, a0) {
  b2 <- a2 / a3; b1 <- a1 / a3; b0 <- a0 / a3
  p <- b1 - b2^2 / 3; q <- 2 * b2^3 / 27 - b2 * b1 / 3 + b0
  disc <- -4 * p^3 - 27 * q^2
  R <- matrix(NA_real_, length(a1), 3)
  three <- disc > 0
  if (any(three)) {
    pp <- p[three]; qq <- q[three]
    ac <- acos(pmin(pmax(3 * qq / (2 * pp) * sqrt(-3 / pp), -1), 1))
    for (k in 0:2) R[three, k + 1] <- 2 * sqrt(-pp / 3) * cos(ac / 3 - 2 * pi * k / 3) - b2 / 3
  }
  one <- !three
  if (any(one)) {
    pp <- p[one]; qq <- q[one]; dt <- sqrt(pmax(qq^2 / 4 + pp^3 / 27, 0))
    u <- -qq / 2 + dt; v <- -qq / 2 - dt
    R[one, 1] <- sign(u) * abs(u)^(1 / 3) + sign(v) * abs(v)^(1 / 3) - b2 / 3
  }
  R
}

# does each row's set {g <= 0} contain c_true? (g = W - threshold on the grid)
contains <- function(G, cgrid, c_true) {
  k <- findInterval(c_true, cgrid)            # c_true in [cgrid[k], cgrid[k+1])
  w <- (c_true - cgrid[k]) / (cgrid[k + 1] - cgrid[k])
  g_at <- (1 - w) * G[, k] + w * G[, k + 1]   # linear interpolation of W - thr
  g_at <= 0
}

set.seed(SEED)
rows <- list()
ops_run <- if (QUICK) ints[[1]]$operator[c(1, 15)] else ints[[1]]$operator
for (op in ops_run) {
  l <- sapply(raw$bins, function(b) raw$param[[b]][[op]])
  q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
  chat <- ints[[1]]$chat[ints[[1]]$operator == op]
  w2 <- ints[[2]]$width_exact[ints[[2]]$operator == op]
  lo <- min(0, chat) - 2.5 * w2; hi <- max(0, chat) + 2.5 * w2
  cgrid <- sort(unique(c(seq(lo, hi, length.out = N_C), 0, chat)))
  # exact thresholds: 150 nodes, spline in between
  cnode <- seq(lo, hi, length.out = if (QUICK) 12 else 150)
  thr_exact <- lapply(LEVELS, function(k) {
    a <- 2 * pnorm(-k)
    v <- vapply(cnode, function(c0) { g <- geometry_at(c0, l, q)
      if (QUICK) quantile_W(g["xi"], g["phi"], a) else
      quantile_W(g["xi"], g["phi"], a, r_max = REF$r_max, n_theta = REF$n_theta, n_r = REF$n_r) }, 0)
    splinefun(cnode, v, method = "natural")
  })
  # spline error at 10 random points
  chk <- runif(if (QUICK) 1 else 10, lo, hi)
  spl_err <- max(abs(vapply(chk, function(c0) { g <- geometry_at(c0, l, q)
    quantile_W(g["xi"], g["phi"], 2 * pnorm(-1), r_max = REF$r_max, n_theta = REF$n_theta, n_r = REF$n_r) }, 0) -
    thr_exact[[1]](chk)))

  dl <- as.numeric(D %*% l); dq <- as.numeric(D %*% q)
  k2c <- sum(l * dl); k3 <- 2 * sum(l * dq); k4 <- sum(q * dq)
  for (c_true in c(0, chat)) {
    mu <- 1 + c_true * l + c_true^2 * q
    Z <- matrix(rnorm(B * length(l)), B) %*% t(Lc)
    R0 <- sweep(Z, 2, mu - 1, "+")                       # Y - 1, one row per dataset
    k0 <- rowSums((R0 %*% D) * R0); k1 <- -2 * as.numeric(R0 %*% dl); k2 <- k2c - 2 * as.numeric(R0 %*% dq)
    dfun <- function(cc) k0 + k1 * cc + k2 * cc^2 + k3 * cc^3 + k4 * cc^4
    Rts <- cubic_roots_vec(4 * k4, 3 * k3, 2 * k2, k1)  # d'(c) = 0
    dmin <- pmin(dfun(Rts[, 1]), dfun(Rts[, 2]), dfun(Rts[, 3]), na.rm = TRUE)
    Wc <- outer(k0, rep(1, length(cgrid))) + outer(k1, cgrid) + outer(k2, cgrid^2) +
          matrix(k3 * cgrid^3 + k4 * cgrid^4, B, length(cgrid), byrow = TRUE) - dmin
    W_true <- dfun(c_true) - dmin
    edge <- mean(Wc[, 1] < 9 | Wc[, length(cgrid)] < 9)  # acceptance set touching the grid ends at any level
    for (k in LEVELS) {
      cv <- 2 * pnorm(k) - 1
      te <- thr_exact[[k]](cgrid); t_true <- thr_exact[[k]](c_true)
      cov_dir <- c(exact = mean(contains(sweep(Wc, 2, te), cgrid, c_true)),
                   chi2 = mean(contains(Wc - qchisq(cv, 1), cgrid, c_true)),
                   chibar = mean(contains(Wc - qchisq(2 * cv - 1, 1), cgrid, c_true)))
      cov_inv <- c(exact = mean(W_true <= t_true), chi2 = mean(W_true <= qchisq(cv, 1)),
                   chibar = mean(W_true <= qchisq(2 * cv - 1, 1)))
      acc <- sweep(Wc, 2, te) <= 0                       # datasets x grid
      ncomp <- colSums(diff(t(acc)) == 1) + acc[, 1]   # number of accepted runs per dataset
      nsplit <- function(A) mean((colSums(diff(t(A)) == 1) + A[, 1]) > 1)
      multi_chi2 <- nsplit(Wc <= qchisq(cv, 1)); multi_chibar <- nsplit(Wc <= qchisq(2 * cv - 1, 1))
      dc <- diff(cgrid)
      len_of <- function(A) as.numeric(((A[, -1] + A[, -ncol(A)]) / 2) %*% dc)
      ctr_of <- function(A) as.numeric(A %*% cgrid) / pmax(rowSums(A), 1)
      sets <- list(exact = acc, chi2 = Wc <= qchisq(cv, 1), chibar = Wc <= qchisq(2 * cv - 1, 1))
      lens <- lapply(sets, len_of); ctrs <- lapply(sets, ctr_of)
      rows[[length(rows) + 1]] <- data.frame(operator = op, c_true = c_true, level_sigma = k,
        nominal = cv, B = B, cover_exact = cov_dir["exact"], cover_chi2 = cov_dir["chi2"],
        cover_chibar = cov_dir["chibar"], inv_exact = cov_inv["exact"], inv_chi2 = cov_inv["chi2"],
        inv_chibar = cov_inv["chibar"], frac_multi_component = mean(ncomp > 1),
        frac_multi_chi2 = multi_chi2, frac_multi_chibar = multi_chibar, frac_edge = edge,
        spline_err = spl_err,
        mean_len_exact = mean(lens$exact), mean_len_chi2 = mean(lens$chi2), mean_len_chibar = mean(lens$chibar),
        med_len_exact = median(lens$exact), med_len_chi2 = median(lens$chi2), med_len_chibar = median(lens$chibar),
        shift_exact = mean(ctrs$exact - c_true), shift_chi2 = mean(ctrs$chi2 - c_true),
        shift_chibar = mean(ctrs$chibar - c_true), row.names = NULL)
      r <- rows[[length(rows)]]
      cat(sprintf("%-6s c_true=%+.3f %dsig: exact %.4f (inv %.4f)  chi2 %.4f  chibar %.4f  [nominal %.4f; multi %.3f; edge %.4f; spline %.1e]\n",
                  op, c_true, k, r$cover_exact, r$inv_exact, r$cover_chi2, r$cover_chibar, cv,
                  r$frac_multi_component, edge, spl_err))
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(root, "tables", if (QUICK) "coverage_length_quick.csv" else "coverage_length.csv"), row.names = FALSE)
cat(sprintf("MC standard error: %.2f pp at 1 sigma, %.2f pp at 2 sigma (B = %d)\n",
            100 * sqrt(0.6827 * 0.3173 / B), 100 * sqrt(0.9545 * 0.0455 / B), B))
