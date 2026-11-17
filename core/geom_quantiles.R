# Tail probabilities and quantiles of W(c0) by quadrature of the bivariate
# standard normal density over {W(c0) > w} in polar coordinates
# (Section 3.4). The standard grid has 48 angular and 400 radial nodes up to
# radius 10; the reference grid uses n_theta = 192, n_r = 800, r_max = 12.
# Requires geom_core.R to be sourced first.

# P(W(c0) > w | xi, phi)
tail_prob <- function(w, xi, phi, r_max = 10.0, n_theta = 48, n_r = 400) {
  thetas <- seq(0, 2 * pi, length.out = n_theta + 1)[1:n_theta]
  r <- seq(1e-6, r_max, length.out = n_r)
  dr <- r[2] - r[1]
  weight <- r * exp(-r^2 / 2)

  zu_mat <- outer(r, cos(thetas))  # n_r x n_theta
  zv_mat <- outer(r, sin(thetas))
  Wv <- compute_w_vectorized(as.vector(zu_mat), as.vector(zv_mat), xi, phi)
  dim(Wv) <- c(n_r, n_theta)

  indicator <- (Wv > w) * weight
  # trapezoidal rule in r for each theta, then average over theta
  col_integrals <- colSums((indicator[-1, , drop = FALSE] + indicator[-n_r, , drop = FALSE]) / 2) * dr
  mean(col_integrals)
}

# upper-alpha quantile of W(c0); 0 if the probability of W > 0 is below alpha
quantile_W <- function(xi, phi, alpha, w_hi0 = 60.0, ...) {
  f0 <- tail_prob(0, xi, phi, ...) - alpha
  if (f0 <= 0) return(0.0)
  w_hi <- w_hi0
  f <- function(w) tail_prob(w, xi, phi, ...) - alpha
  fhi <- f(w_hi)
  tries <- 0
  while (fhi > 0 && tries < 6) {
    w_hi <- w_hi * 2
    fhi <- f(w_hi)
    tries <- tries + 1
  }
  uniroot(f, c(0, w_hi), tol = 1e-4)$root
}
