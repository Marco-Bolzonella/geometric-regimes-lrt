# Poisson simulation for Table S2 and Figure 6. Counts are
# independent Poisson with means m (1 + c l_i + c^2 q_i). For each dataset
# the statistic for c = 0 is computed from the Gaussian working likelihood
# with covariance I/m for Y/m, the statistic the method calibrates; the
# minimum over c is found on a grid and refined by golden-section search.

# golden-section minimisation of f (one value per row) on brackets [lo, hi]
golden_vec <- function(f, lo, hi, iters = 60) {
  g <- (sqrt(5) - 1) / 2
  a <- lo; b <- hi
  x1 <- b - g * (b - a); x2 <- a + g * (b - a)
  f1 <- f(x1); f2 <- f(x2)
  for (it in seq_len(iters)) {
    left <- f1 < f2
    b <- ifelse(left, x2, b); a <- ifelse(left, a, x1)
    x2n <- ifelse(left, x1, a + g * (b - a)); x1n <- ifelse(left, b - g * (b - a), x2)
    x1 <- x1n; x2 <- x2n
    f1 <- f(x1); f2 <- f(x2)
  }
  (a + b) / 2
}

# B datasets, processed in chunks to limit memory
simulate_poisson_W <- function(l, q, mu_scale, B, c_grid, seed, chunk_size = 4000) {
  set.seed(seed)
  n_bins <- length(l)
  n_grid <- length(c_grid)

  Rel <- outer(rep(1, n_bins), rep(1, n_grid)) + outer(l, c_grid) + outer(q, c_grid^2)
  if (any(Rel <= 0)) {
    stop("mu_i(c) <= 0 somewhere on the grid -- narrow c_grid or check l,q")
  }

  # Gaussian working likelihood: d(c) = m || Y/m - (1 + c l + c^2 q) ||^2
  colsum_Rel2 <- colSums(Rel^2)
  n_chunks <- ceiling(B / chunk_size)
  W_all <- numeric(B)
  boundary_hits <- 0

  for (ch in seq_len(n_chunks)) {
    b_lo <- (ch - 1) * chunk_size + 1
    b_hi <- min(ch * chunk_size, B)
    n_b <- b_hi - b_lo + 1

    Y <- matrix(rpois(n_b * n_bins, lambda = mu_scale), nrow = n_b, ncol = n_bins)
    Yr <- Y / mu_scale
    dG <- mu_scale * (rowSums(Yr^2) - 2 * (Yr %*% Rel) + matrix(colsum_Rel2, n_b, n_grid, byrow = TRUE))
    gi <- max.col(-dG, ties.method = "first")
    boundary_hits <- boundary_hits + sum(gi %in% c(1, n_grid))

    # refine on [c_{k-1}, c_{k+1}] around the best grid point
    lo_g <- c_grid[pmax(gi - 1, 1)]; hi_g <- c_grid[pmin(gi + 1, n_grid)]
    dfun <- function(cc) mu_scale * rowSums((Yr - (1 + outer(cc, l) + outer(cc^2, q)))^2)
    cg <- golden_vec(dfun, lo_g, hi_g)
    dmin <- pmin(dfun(cg), dG[cbind(seq_len(n_b), gi)])
    d0 <- mu_scale * rowSums((Yr - 1)^2)
    W_all[b_lo:b_hi] <- d0 - dmin
  }

  boundary_frac <- boundary_hits / B
  if (boundary_frac > 0.001) {
    warning(sprintf("simulate_poisson_W: MLE hit the c_grid boundary for %.2f%% of replicates -- widen c_grid",
                     boundary_frac * 100))
  }
  list(W = W_all, boundary_frac = boundary_frac)
}
