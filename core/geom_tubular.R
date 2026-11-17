# Density of W(c0) from the tubular-coordinate integral, Eq. (15):
#
#   f_W(w) = int f_Z(Phi(tau, s)) |J(tau, s)| / (2 |h(tau)|) dtau,
#   Phi(tau, s) = gamma(tau) + s n(tau),  s = (w - |gamma(tau)|^2) / (2 h(tau)),
#
# with h, kappa and J as in Lemma C.1, restricted to points at which tau is
# the global minimiser. Used only to check geom_quantiles.R; the minimiser
# test uses compute_w_general() rather than the vectorised solver. The
# integral in tau is split at tau = 0, at the cusp -cos(phi)/xi, at the
# roots of |gamma(tau)|^2 = w and wherever the minimiser test changes.

tub_parts <- function(tau, xi, phi) {
  u <- tau + 0.5 * xi * cos(phi) * tau^2
  v <- 0.5 * xi * sin(phi) * tau^2
  du <- 1 + xi * cos(phi) * tau
  dv <- xi * sin(phi) * tau
  sp <- sqrt(du^2 + dv^2)
  nu <- -dv / sp; nv <- du / sp   # unit normal
  h <- u * nu + v * nv
  kappa <- xi * sin(phi) / sp^3
  list(u = u, v = v, nu = nu, nv = nv, sp = sp, h = h, kappa = kappa, g2 = u^2 + v^2)
}

# is tau the global minimiser of the distance from Phi(tau, s)?
is_global <- function(tau, s, xi, phi, P) {
  zu <- P$u + s * P$nu; zv <- P$v + s * P$nv
  taus <- cubic_roots(zu, zv, xi, phi)
  if (length(taus) <= 1) return(TRUE)
  cand <- c(taus[1], taus[length(taus)])
  d2 <- d2_dist(cand, zu, zv, xi, phi)
  # tau is global iff it attains the minimum distance s^2
  min(d2) >= s^2 - 1e-10 * max(1, s^2)
}

tub_integrand1 <- function(tau, w, xi, phi) {
  P <- tub_parts(tau, xi, phi)
  if (abs(P$h) < 1e-300) return(0)
  s <- (w - P$g2) / (2 * P$h)
  if (!is_global(tau, s, xi, phi, P)) return(0)
  zu <- P$u + s * P$nu; zv <- P$v + s * P$nv
  J <- P$sp * (1 - s * P$kappa)
  exp(-(zu^2 + zv^2) / 2) / (2 * pi) * abs(J) / (2 * abs(P$h))
}
tub_integrand <- function(tau, w, xi, phi) vapply(tau, tub_integrand1, numeric(1), w = w, xi = xi, phi = phi)

# f_W(w) by adaptive quadrature in tau
fW_tubular <- function(w, xi, phi, tau_max = NULL, n_scan = 4001, rel.tol = 1e-10) {
  # beyond |gamma(tau)| = sqrt(w) + 12 the integrand is negligible
  if (is.null(tau_max)) {
    R <- sqrt(w) + 12
    tau_max <- max(R, sqrt(2 * R / (xi * max(sin(phi), 1e-3))) + 1)
  }
  # roots of |gamma(tau)|^2 = w: quartic (tau + a tau^2)^2 + (b tau^2)^2 = w
  a <- 0.5 * xi * cos(phi); b <- 0.5 * xi * sin(phi)
  pr <- polyroot(c(-w, 0, 1, 2 * a, a^2 + b^2))
  pr <- Re(pr[abs(Im(pr)) < 1e-8])
  # scan grid, refined near tau = 0 and near the roots above
  lg <- 10^seq(-9, log10(tau_max), length.out = 400)
  near <- unlist(lapply(pr, function(r) r + c(-1, 1) %o% 10^seq(-9, 0, length.out = 120)))
  grid <- c(seq(-tau_max, tau_max, length.out = n_scan), 0, lg, -lg, near)
  grid <- sort(unique(grid[abs(grid) <= tau_max]))
  brk <- c(0, -cos(phi) / xi, pr[abs(pr) < tau_max])
  # indicator switches
  ind <- vapply(grid, function(t) tub_integrand1(t, w, xi, phi) > 0, logical(1))
  sw <- which(diff(ind) != 0)
  for (k in sw) {
    lo <- grid[k]; hi <- grid[k + 1]; il <- ind[k]
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      if ((tub_integrand1(mid, w, xi, phi) > 0) == il) lo <- mid else hi <- mid
    }
    brk <- c(brk, (lo + hi) / 2)
  }
  brk <- sort(unique(c(-tau_max, brk[abs(brk) < tau_max], tau_max)))
  total <- 0
  for (k in seq_len(length(brk) - 1)) {
    a1 <- brk[k]; b1 <- brk[k + 1]
    if (b1 - a1 < 1e-14) next
    r <- integrate(tub_integrand, a1, b1, w = w, xi = xi, phi = phi,
                   rel.tol = rel.tol, abs.tol = 0, subdivisions = 2000L, stop.on.error = FALSE)
    total <- total + r$value
  }
  total
}

gauss_legendre <- function(n) {
  # Golub--Welsch: nodes and weights on [-1, 1]
  k <- seq_len(n - 1); b <- k / sqrt(4 * k^2 - 1)
  e <- eigen(diag(0, n) + rbind(cbind(0, diag(b, n - 1)), 0) + t(rbind(cbind(0, diag(b, n - 1)), 0)),
             symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# P(W <= w): Gauss-Legendre in x = sqrt(t) on (0, sqrt(w)]
cdf_tubular <- function(w, xi, phi, n_gl = 20, n_pieces = 4) {
  gl <- gauss_legendre(n_gl)
  edges <- seq(0, sqrt(w), length.out = n_pieces + 1)
  total <- 0
  for (k in seq_len(n_pieces)) {
    a <- edges[k]; b <- edges[k + 1]
    x <- (b - a) / 2 * gl$x + (a + b) / 2
    fx <- vapply(x, function(xx) fW_tubular(xx^2, xi, phi) * 2 * xx, numeric(1))
    total <- total + (b - a) / 2 * sum(gl$w * fx)
  }
  total
}

# P(W > w): Gauss-Legendre in x = sqrt(t) on [w, w + 80]. Preferred to
# cdf_tubular() when xi is large, because W then puts about half its mass
# within O(1/xi) of zero.
tail_tubular <- function(w, xi, phi, n_gl = 20, n_pieces = 8, span = 80) {
  gl <- gauss_legendre(n_gl)
  edges <- seq(sqrt(w), sqrt(w + span), length.out = n_pieces + 1)
  total <- 0
  for (k in seq_len(n_pieces)) {
    a <- edges[k]; b <- edges[k + 1]
    x <- (b - a) / 2 * gl$x + (a + b) / 2
    fx <- vapply(x, function(xx) fW_tubular(xx^2, xi, phi) * 2 * xx, numeric(1))
    total <- total + (b - a) / 2 * sum(gl$w * fx)
  }
  total
}
