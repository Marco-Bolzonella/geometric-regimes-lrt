# Figure 8 (Section 5.3): (xi(c0), phi(c0)) across the 68.3% acceptance set
# on the measured data (needs tables/measured_intervals_1.csv). The maximum
# of xi(c0) is at c* = -<T(0), N> / ||N||^2, where T(c0) is shortest. Also
# prints the ranges of xi(c0) quoted in the text and the partial-correction
# thresholds for cQd1.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "plot_theme.R"))

raw <- load_raw(file.path(root, "data"))
M <- whitening_matrix(raw$cov)
ends <- read.csv(file.path(root, "tables", "measured_intervals_1.csv"))
OPS <- c("ctj8", "cQu8", "cG", "cQd1")
ALPHA1 <- 2 * pnorm(-1)

vecs <- function(op) {
  l <- sapply(raw$bins, function(b) raw$param[[b]][[op]])
  q <- sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]])
  list(T0 = as.numeric(M %*% l), N = as.numeric(2 * (M %*% q)))
}
geom_at <- function(op, c0) {
  v <- vecs(op)
  t(vapply(c0, function(c) {
    Tc <- v$T0 + c * v$N; e1 <- Tc / sqrt(sum(Tc^2)); a <- sum(v$N * e1)
    c(xi = sqrt(sum(v$N^2)) / sum(Tc^2),
      phi_deg = atan2(sqrt(sum((v$N - a * e1)^2)), a) * 180 / pi)
  }, numeric(2)))
}
c_star <- function(op) { v <- vecs(op); -sum(v$T0 * v$N) / sum(v$N^2) }

# "[a,b]U[c,d]" -> 2-column matrix
parse_set <- function(s) {
  m <- regmatches(s, gregexpr("-?[0-9.]+", s))[[1]]
  matrix(as.numeric(m), ncol = 2, byrow = TRUE)
}

curve <- list(); pts <- list()
for (op in OPS) {
  e <- ends[ends$operator == op, ]
  set <- parse_set(e$set_exact)
  lo <- min(set); hi <- max(set); cs <- c_star(op)
  # denser sampling near c*
  c0 <- seq(lo, hi, length.out = 2000)
  if (cs > lo && cs < hi) {
    t <- seq(-1, 1, length.out = 3000)^3
    c0 <- c(c0, cs + ifelse(t < 0, -t * (lo - cs), t * (hi - cs)), cs)
  }
  c0 <- sort(unique(c0))
  g <- geom_at(op, c0)
  inside <- vapply(c0, function(c) any(c >= set[, 1] - 1e-9 & c <= set[, 2] + 1e-9), logical(1))
  # excluded stretches are drawn dashed
  seg <- cumsum(c(TRUE, diff(inside) != 0))
  curve[[op]] <- data.frame(operator = op, c0 = c0, g, accepted = inside, seg = paste(op, seg))

  # range of xi(c0) on each component
  for (r in seq_len(nrow(set))) {
    a <- set[r, 1]; b <- set[r, 2]
    cand <- c(a, b, if (cs > a && cs < b) cs)
    xs <- geom_at(op, cand)[, "xi"]
    cc <- seq(a, b, length.out = 200001)
    frac <- mean(geom_at(op, cc[seq(1, length(cc), by = 20)])[, "xi"] > 100)
    cat(sprintf("%-5s component %d [%.4f, %.4f]: xi from %.2f to %.2f; xi > 100 on %.1f%%\n",
                op, r, a, b, min(xs), max(xs), 100 * frac))
  }
  mk <- data.frame(operator = op,
                   c0 = c(e$chat, e$lo_chibar, e$hi_chibar, as.vector(t(set))),
                   what = c("MLE", "chibar2 endpoint", "chibar2 endpoint",
                            rep("Exact endpoint", length(set))))
  pts[[op]] <- cbind(mk, geom_at(op, mk$c0))
}
cdf <- do.call(rbind, curve); pdf_ <- do.call(rbind, pts)
cdf$operator  <- factor(cdf$operator,  levels = OPS)
pdf_$operator <- factor(pdf_$operator, levels = OPS)
pdf_$what <- factor(pdf_$what, levels = c("MLE", "chibar2 endpoint", "Exact endpoint"))
print(pdf_, digits = 4)

# partial corrections at the chibar2 endpoints of cQd1
e <- ends[ends$operator == "cQd1", ]
g0 <- geom_at("cQd1", 0)
pc_rows <- list()
for (cb in c(e$lo_chibar, e$hi_chibar)) {
  g <- geom_at("cQd1", cb)
  full <- quantile_W(g[1, "xi"], g[1, "phi_deg"] * pi / 180, ALPHA1, r_max = 12, n_theta = 192, n_r = 800)
  xi_only <- quantile_W(g[1, "xi"], g0[1, "phi_deg"] * pi / 180, ALPHA1, r_max = 12, n_theta = 192, n_r = 800)
  phi_only <- quantile_W(g0[1, "xi"], g[1, "phi_deg"] * pi / 180, ALPHA1, r_max = 12, n_theta = 192, n_r = 800)
  cat(sprintf("cQd1 chibar2 endpoint c0=%+.4f: xi=%.3f phi=%.1f  exact thr=%.3f  xi adapted only=%.3f  phi adapted only=%.3f  (chibar2 thr=%.3f)\n",
              cb, g[1, "xi"], g[1, "phi_deg"], full, xi_only, phi_only, qchisq(1 - 2 * ALPHA1, 1)))
  pc_rows[[length(pc_rows) + 1]] <- data.frame(operator = "cQd1", c0 = cb, xi = g[1, "xi"],
    phi_deg = g[1, "phi_deg"], thr_exact = full, thr_xi_adapted_only = xi_only,
    thr_phi_adapted_only = phi_only, thr_chibar = qchisq(1 - 2 * ALPHA1, 1))
}
write.csv(do.call(rbind, pc_rows), file.path(root, "tables", "partial_correction.csv"), row.names = FALSE)

what_lab <- c("MLE" = expression(hat(c)),
              "chibar2 endpoint" = expression(bar(chi)^2 ~ "endpoint"),
              "Exact endpoint" = "Exact endpoint")
p <- ggplot(cdf, aes(xi, phi_deg)) +
  geom_path(aes(group = seg, linetype = accepted), colour = COL_CURVE, linewidth = 0.5,
            show.legend = FALSE) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22")) +
  geom_point(data = pdf_, aes(shape = what, fill = what), size = 1.8, stroke = 0.4) +
  scale_shape_manual(values = c(21, 24, 22), labels = what_lab, name = NULL) +
  scale_fill_manual(values = c(COL_EXACT, COL_CHIBAR, "white"), labels = what_lab, name = NULL) +
  scale_x_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, big.mark = ",")) +
  scale_y_continuous(breaks = seq(0, 180, 45), limits = c(0, 180)) +
  facet_wrap(~operator, nrow = 1, scales = "free_x") +
  labs(x = expression(xi(c[0]) ~ "(log scale)"), y = expression(varphi(c[0]) ~ "(degrees)")) +
  theme_paper() + theme(legend.position = "bottom")

save_fig(p, file.path(root, "figures", "trajectory_plot.png"), height = 2.9)
