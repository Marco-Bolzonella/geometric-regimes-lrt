# Likelihood ratio statistic W(c0) for a realisation (zu, zv) of Z_L and
# geometry (xi, phi), from the real roots of the cubic stationarity
# equation (11). compute_w_general() uses polyroot() and serves as the
# reference; compute_w_vectorized() solves the cubic in closed form and is
# used everywhere else. self_test() compares the two.

# real roots of (1/2) xi^2 tau^3 + (3/2) xi cos(phi) tau^2
#   + (1 - xi cos(phi) zu - xi sin(phi) zv) tau - zu = 0
cubic_roots <- function(zu, zv, xi, phi) {
  A <- 0.5 * xi^2
  B <- 1.5 * xi * cos(phi)
  C <- 1 - xi * cos(phi) * zu - xi * sin(phi) * zv
  D <- -zu
  if (A == 0) {
    coeffs <- if (C != 0) c(D, C) else c(D, 1e-12)
    roots <- polyroot(coeffs)
  } else {
    roots <- polyroot(c(D, C, B, A))
  }
  re <- Re(roots)
  im <- Im(roots)
  sort(re[abs(im) < 1e-7 * pmax(1, abs(re))])
}

# squared distance from z to the curve point at tau
d2_dist <- function(tau, zu, zv, xi, phi) {
  u <- tau + 0.5 * xi * cos(phi) * tau^2
  v <- 0.5 * xi * sin(phi) * tau^2
  (zu - u)^2 + (zv - v)^2
}

# W(c0) for a single (zu, zv); with three real roots the minimiser is one
# of the two outer ones (Lemma B.2)
compute_w_general <- function(zu, zv, xi, phi) {
  taus <- cubic_roots(zu, zv, xi, phi)
  if (length(taus) == 0) {
    grid <- seq(-50, 50, length.out = 20001)
    d2v <- d2_dist(grid, zu, zv, xi, phi)
    tau_star <- grid[which.min(d2v)]
  } else if (length(taus) == 1) {
    tau_star <- taus[1]
  } else {
    cand <- c(taus[1], taus[length(taus)])
    d2v <- d2_dist(cand, zu, zv, xi, phi)
    tau_star <- cand[which.min(d2v)]
  }
  zu^2 + zv^2 - d2_dist(tau_star, zu, zv, xi, phi)
}

# W(c0) for vectors zu, zv: trigonometric solution when the cubic has
# three real roots, Cardano's formula otherwise
compute_w_vectorized <- function(zu, zv, xi, phi) {
  n <- length(zu)
  A <- 0.5 * xi^2
  B <- 1.5 * xi * cos(phi)

  if (A < 1e-12) {
    tau0 <- zu
    u <- tau0
    v <- rep(0, n)
    return(zu^2 + zv^2 - ((zu - u)^2 + (zv - v)^2))
  }

  C <- 1 - xi * cos(phi) * zu - xi * sin(phi) * zv
  D <- -zu

  a2 <- B / A
  a1 <- C / A
  a0 <- D / A

  p <- a1 - a2^2 / 3
  q <- (2 * a2^3) / 27 - (a2 * a1) / 3 + a0
  disc <- -4 * p^3 - 27 * q^2

  tau_candidates <- matrix(NA_real_, nrow = 3, ncol = n)

  three_root <- disc > 0
  if (any(three_root)) {
    pp <- p[three_root]
    qq <- q[three_root]
    arg <- pmin(pmax((3 * qq) / (2 * pp) * sqrt(-3 / pp), -1), 1)
    acos_term <- acos(arg)
    for (k in 0:2) {
      y_k <- 2 * sqrt(-pp / 3) * cos(acos_term / 3 - 2 * pi * k / 3)
      tau_candidates[k + 1, three_root] <- y_k - a2 / 3
    }
  }

  one_root <- !three_root
  if (any(one_root)) {
    pp <- p[one_root]
    qq <- q[one_root]
    disc_term <- sqrt(pmax(qq^2 / 4 + pp^3 / 27, 0))
    inner1 <- -qq / 2 + disc_term
    inner2 <- -qq / 2 - disc_term
    y <- sign(inner1) * abs(inner1)^(1 / 3) + sign(inner2) * abs(inner2)^(1 / 3)
    tau_candidates[1, one_root] <- y - a2 / 3
  }

  d2_of <- function(tau, zu_, zv_) {
    u <- tau + 0.5 * xi * cos(phi) * tau^2
    v <- 0.5 * xi * sin(phi) * tau^2
    (zu_ - u)^2 + (zv_ - v)^2
  }

  # outer roots (equal when there is only one)
  tau_lo <- pmin(tau_candidates[1, ], tau_candidates[2, ], tau_candidates[3, ], na.rm = TRUE)
  tau_hi <- pmax(tau_candidates[1, ], tau_candidates[2, ], tau_candidates[3, ], na.rm = TRUE)

  d2_lo <- d2_of(tau_lo, zu, zv)
  d2_hi <- d2_of(tau_hi, zu, zv)
  zu^2 + zv^2 - pmin(d2_lo, d2_hi)
}

# W(c0) at polar coordinates (r, theta)
compute_w_radial <- function(r, theta, xi, phi) {
  zu <- r * cos(theta)
  zv <- r * sin(theta)
  compute_w_vectorized(zu, zv, xi, phi)
}

# compare the two implementations at random points
self_test <- function(seed = 20260909, n = 300, verbose = TRUE) {
  set.seed(seed)
  max_err <- 0
  for (i in seq_len(n)) {
    zu <- rnorm(1) * 3
    zv <- rnorm(1) * 3
    xi <- runif(1, 0.05, 3000)
    phi <- runif(1, 0.01, pi - 0.01)
    w_ref <- compute_w_general(zu, zv, xi, phi)
    w_fast <- compute_w_vectorized(zu, zv, xi, phi)
    max_err <- max(max_err, abs(w_ref - w_fast))
  }
  if (verbose) {
    cat(sprintf("self-test: max |compute_w_general - compute_w_vectorized| over %d random points = %.3e\n",
                n, max_err))
  }
  stopifnot(max_err < 1e-8)
  invisible(max_err)
}

if (sys.nframe() == 0) {
  self_test()
}
