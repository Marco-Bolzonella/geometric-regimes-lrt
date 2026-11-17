# Table S1, last column: coverage at the reference-grid quantile
# recomputed from Eq. (15), for two (xi, phi) per regime at 1-4 sigma.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "geom_tubular.R"))
source(file.path(root, "core", "load_operators.R"))

ops <- build_table1(file.path(root, "data"))
rownames(ops) <- ops$operator
POINTS <- c("ctGRe", "cQu8", "cQj38", "cQj31", "cG", "cQd1")

tail_ref <- function(w, xi, phi) tail_prob(w, xi, phi, r_max = 12, n_theta = 192, n_r = 800)

rows <- list()
for (op in POINTS) {
  xi <- ops[op, "xi"]; phi <- ops[op, "phi_rad"]
  for (k in 1:4) {
    cover <- pnorm(k) - pnorm(-k)
    f <- function(w) tail_ref(w, xi, phi) - (1 - cover)
    w_k <- uniroot(f, c(1e-6, 40), tol = 1e-10)$root
    t0 <- Sys.time()
    F_tub <- 1 - tail_tubular(w_k, xi, phi)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    F_std <- 1 - tail_prob(w_k, xi, phi)
    F_ref <- 1 - tail_ref(w_k, xi, phi)
    rows[[length(rows) + 1]] <- data.frame(
      operator = op, xi = xi, phi_deg = phi * 180 / pi, k_sigma = k, nominal = cover,
      w_k = w_k, cover_tubular = F_tub, cover_plane_std = F_std, cover_plane_ref = F_ref,
      diff_std_pp = 100 * (F_std - F_tub), diff_ref_pp = 100 * (F_ref - F_tub),
      seconds_tubular = secs)
    cat(sprintf("%-6s xi=%9.2f phi=%6.1f %ds  w=%.4f  tub=%.6f  std=%.6f (%+.4f pp)  ref=%.6f (%+.4f pp)  [%.0fs]\n",
                op, xi, phi * 180 / pi, k, w_k, F_tub, F_std, 100 * (F_std - F_tub),
                F_ref, 100 * (F_ref - F_tub), secs))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(root, "tables", "tubular_check.csv"), row.names = FALSE)
cat(sprintf("max |plane std - tubular| = %.4f pp; max |plane ref - tubular| = %.4f pp\n",
            max(abs(res$diff_std_pp)), max(abs(res$diff_ref_pp))))
