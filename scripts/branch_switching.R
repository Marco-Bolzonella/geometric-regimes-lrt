# Branch switching at xi = 0.53 (Section 4.2 and Appendix F): splits the 1-sigma coverage
# error of chi2_1 into a local part, computed with the local minimiser of the
# squared distance closest to z_u (the minimiser as xi -> 0), and the
# remainder, due to realisations whose global minimiser is the other local
# minimiser. The middle root of the cubic, a local maximum, is excluded.
# Polar quadrature, 720 x 3000 nodes.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))

th <- (seq_len(720) - 0.5) * 2 * pi / 720
r <- (seq_len(3000) - 0.5) * 12 / 3000
g <- expand.grid(r = r, th = th)
zu <- g$r * cos(g$th); zv <- g$r * sin(g$th)
wt <- g$r * exp(-g$r^2 / 2) * (12 / 3000) * (2 * pi / 720) / (2 * pi)
cov1 <- 2 * pnorm(1) - 1; w1 <- qchisq(cov1, 1)

# real roots of the cubic, NA where there is only one
real_roots <- function(zu, zv, xi, phi) {
  A <- 0.5 * xi^2; B <- 1.5 * xi * cos(phi)
  C <- 1 - xi * cos(phi) * zu - xi * sin(phi) * zv; D <- -zu
  a2 <- B / A; a1 <- C / A; a0 <- D / A
  p <- a1 - a2^2 / 3; q <- 2 * a2^3 / 27 - a2 * a1 / 3 + a0
  three <- (-4 * p^3 - 27 * q^2) > 0
  R <- matrix(NA_real_, length(zu), 3)
  if (any(three)) {
    pp <- p[three]; qq <- q[three]
    ac <- acos(pmin(pmax(3 * qq / (2 * pp) * sqrt(-3 / pp), -1), 1))
    for (k in 0:2) R[three, k + 1] <- 2 * sqrt(-pp / 3) * cos(ac / 3 - 2 * pi * k / 3) - a2 / 3
  }
  one <- !three
  pp <- p[one]; qq <- q[one]; dt <- sqrt(pmax(qq^2 / 4 + pp^3 / 27, 0))
  i1 <- -qq / 2 + dt; i2 <- -qq / 2 - dt
  R[one, 1] <- sign(i1) * abs(i1)^(1 / 3) + sign(i2) * abs(i2)^(1 / 3) - a2 / 3
  R
}

rows <- list()
for (case in list(c(0.53, 20), c(0.53, 127), c(0.53, 90))) {
  xi <- case[1]; ph <- case[2] * pi / 180
  W <- compute_w_vectorized(zu, zv, xi, ph)
  R <- real_roots(zu, zv, xi, ph)
  # keep local minima only: d^2''(tau) > 0, i.e. the derivative of the cubic is positive
  a <- xi * cos(ph); b <- xi * sin(ph)
  fp <- 1.5 * xi^2 * R^2 + 3 * a * R + (1 - a * zu - b * zv)
  dist <- abs(R - zu); dist[is.na(dist) | fp <= 0] <- Inf
  tl <- R[cbind(seq_len(nrow(R)), max.col(-dist, ties.method = "first"))]
  Wl <- zu^2 + zv^2 - d2_dist(tl, zu, zv, xi, ph)
  tot <- 100 * (sum((W <= w1) * wt) - cov1)
  loc <- 100 * (sum((Wl <= w1) * wt) - cov1)
  psw <- sum((abs(W - Wl) > 1e-9) * wt)
  rows[[length(rows) + 1]] <- data.frame(xi = xi, phi_deg = case[2], total_pp = tot, local_pp = loc,
                                         switching_pp = tot - loc, prob_switch = psw,
                                         recross_height = 2 * sin(ph) / (xi * cos(ph)^2))
  cat(sprintf("xi=%.2f phi=%3d: total %.2f pp = local %.2f + switching %.2f; P(switch) = %.4f\n",
              xi, case[2], tot, loc, tot - loc, psw))
}
write.csv(do.call(rbind, rows), file.path(root, "tables", "branch_switching.csv"), row.names = FALSE)
cat(sprintf("grid mass %.7f\n", sum(wt)))
