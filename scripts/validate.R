# Self-tests: the two W(c0) solvers, and the quadrature against Monte
# Carlo and against the chi2_1 and chibar2 limits.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))

SEED <- 20260909

check_vs_special_case_phi90 <- function() {
  set.seed(SEED)
  max_err <- 0
  for (i in 1:200) {
    zu <- rnorm(1) * 2
    zv <- rnorm(1) * 2
    xi <- runif(1, 0.1, 50)
    phi <- pi / 2
    w_general <- compute_w_general(zu, zv, xi, phi)
    A <- 0.5 * xi^2; C <- 1 - xi * zv; D <- -zu
    roots <- polyroot(c(D, C, 0, A))
    re <- Re(roots); im <- Im(roots)
    real_roots <- re[abs(im) < 1e-7]
    d2v <- (zu - real_roots)^2 + (zv - 0.5 * xi * real_roots^2)^2
    tau_star <- real_roots[which.min(d2v)]
    w_special <- zu^2 + zv^2 - ((zu - tau_star)^2 + (zv - 0.5 * xi * tau_star^2)^2)
    max_err <- max(max_err, abs(w_general - w_special))
  }
  cat(sprintf("[phi=pi/2 special case] max abs error over 200 points: %.2e\n", max_err))
  stopifnot(max_err < 1e-10)
}

check_wilks_limit <- function() {
  alpha <- 0.05
  target <- qchisq(1 - alpha, df = 1)
  for (xi in c(0.01, 0.001)) {
    q <- quantile_W(xi, pi / 2, alpha)
    cat(sprintf("[xi->0 Wilks limit] xi=%s: q=%.4f vs target %.4f (err=%+.4f)\n",
                xi, q, target, q - target))
  }
}

check_chernoff_limit <- function() {
  alpha <- 0.05
  target <- qchisq(1 - 2 * alpha, df = 1)
  for (xi in c(2000, 20000, 200000)) {
    q <- quantile_W(xi, pi / 2, alpha)
    cat(sprintf("[xi->inf Chernoff limit] xi=%s: q=%.4f vs target %.4f (err=%+.4f)\n",
                xi, q, target, q - target))
  }
}

check_vs_montecarlo <- function() {
  set.seed(SEED)
  n_mc <- 200000
  cat("[tail_prob vs Monte Carlo]\n")
  configs <- list(c(30, 90), c(30, 45), c(6, 90), c(150, 60))
  for (cfg in configs) {
    xi <- cfg[1]; phi_deg <- cfg[2]
    phi <- phi_deg * pi / 180
    zu <- rnorm(n_mc); zv <- rnorm(n_mc)
    Wmc <- compute_w_vectorized(zu, zv, xi, phi)
    for (w in c(3.0, 5.0)) {
      p_mc <- mean(Wmc > w)
      se_mc <- sqrt(p_mc * (1 - p_mc) / n_mc)
      p_quad <- tail_prob(w, xi, phi)
      flag <- if (abs(p_quad - p_mc) < 4 * se_mc) "OK" else "CHECK"
      cat(sprintf("  xi=%4g phi=%3gdeg w=%s: MC=%.4f(se=%.4f) quad=%.4f [%s]\n",
                  xi, phi_deg, w, p_mc, se_mc, p_quad, flag))
    }
  }
}

cat("=== 1. compute_w_general vs compute_w_vectorized ===\n")
self_test()
cat("\n=== 2. phi=pi/2 special-case cross-check ===\n")
check_vs_special_case_phi90()
cat("\n=== 3. xi->0 Wilks limit ===\n")
check_wilks_limit()
cat("\n=== 4. xi->inf Chernoff limit ===\n")
check_chernoff_limit()
cat("\n=== 5. tail_prob vs Monte Carlo ===\n")
check_vs_montecarlo()
cat("\nAll checks completed.\n")
