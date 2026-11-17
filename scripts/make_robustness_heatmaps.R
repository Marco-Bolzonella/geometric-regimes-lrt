# Figures 5 and S1 (Section 4.3): error rate at 1 and 2 sigma of the tests based on the three quantiles when
# the coordinates of Z_L are Normal, t5 or t3 with unit variance. Rows of
# the grid are processed in batches and checkpointed; rerun until the
# script prints COMPLETE (run_all.R does this).

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "plot_theme.R"))

SEED <- 20260909
N_XI <- 26
N_PHI <- 22
N_MC <- 400000
LEVELS <- c(1, 2)                                  # k sigma
ALPHAS <- 2 * pnorm(-LEVELS)                       # 0.3173, 0.0455
METHODS <- c("exact", "chi2", "chibar")
CHECKPOINT <- file.path(root, "tables", "robustness_grid_levels_checkpoint.rds")
ROWS_PER_RUN <- 22

set.seed(SEED)

xi_grid <- 10^seq(log10(0.1), log10(20000), length.out = N_XI)
phi_grid <- seq(10, 170, length.out = N_PHI)

scenarios <- c("Normal", "t5", "t3")
sample_noise <- function(scenario, n) {
  if (scenario == "Normal") {
    list(zu = rnorm(n), zv = rnorm(n))
  } else if (scenario == "t5") {
    list(zu = rt(n, 5) * sqrt(3 / 5), zv = rt(n, 5) * sqrt(3 / 5))
  } else if (scenario == "t3") {
    list(zu = rt(n, 3) * sqrt(1 / 3), zv = rt(n, 3) * sqrt(1 / 3))
  }
}

w_chi2_fixed <- qchisq(1 - ALPHAS, df = 1)
w_chibar_fixed <- qchisq(1 - 2 * ALPHAS, df = 1)

if (file.exists(CHECKPOINT)) {
  cp <- readRDS(CHECKPOINT)
  results <- cp$results
  next_row <- cp$next_row
  cat(sprintf("resuming from checkpoint: row %d/%d already done\n", next_row - 1, N_PHI))
} else {
  empty <- matrix(NA_real_, nrow = N_PHI, ncol = N_XI)
  results <- lapply(seq_along(LEVELS), function(k)
    setNames(lapply(scenarios, function(s)
      setNames(lapply(METHODS, function(m) empty), METHODS)), scenarios))
  next_row <- 1
}

last_row <- min(next_row + ROWS_PER_RUN - 1, N_PHI)
if (next_row > N_PHI) {
  cat("COMPLETE -- all rows already done, skipping to plot\n")
} else {
  for (i in next_row:last_row) {
    set.seed(SEED + i)  # per-row seed: results do not depend on how rows are split across runs
    phi <- phi_grid[i] * pi / 180
    for (j in seq_along(xi_grid)) {
      xi <- xi_grid[j]
      w_thr <- vapply(ALPHAS, function(a) quantile_W(xi, phi, a), numeric(1))
      for (s in scenarios) {
        noise <- sample_noise(s, N_MC)
        Wv <- compute_w_vectorized(noise$zu, noise$zv, xi, phi)
        for (k in seq_along(LEVELS)) {
          results[[k]][[s]]$exact[i, j] <- mean(Wv > w_thr[k])
          results[[k]][[s]]$chi2[i, j] <- mean(Wv > w_chi2_fixed[k])
          results[[k]][[s]]$chibar[i, j] <- mean(Wv > w_chibar_fixed[k])
        }
      }
    }
    cat(sprintf("  row %d/%d done\n", i, N_PHI))
  }
  saveRDS(list(results = results, next_row = last_row + 1), CHECKPOINT)
}

if (last_row < N_PHI) {
  cat(sprintf("PARTIAL -- %d/%d rows done, rerun this script to continue\n", last_row, N_PHI))
  quit(status = 0)
}

cat("COMPLETE -- all rows done, producing plot\n")
saveRDS(list(xi_grid = xi_grid, phi_grid = phi_grid, results = results, levels = LEVELS, alphas = ALPHAS),
        file.path(root, "tables", "robustness_grid_levels.rds"))

operators <- build_table1(file.path(root, "data"))

for (k in seq_along(LEVELS)) for (m in METHODS) for (s in scenarios) {
  v <- 100 * results[[k]][[s]][[m]]
  cat(sprintf("%d sigma (nominal %.2f%%) %-6s %-6s: size [%.2f, %.2f] %%\n",
              LEVELS[k], 100 * ALPHAS[k], m, s, min(v), max(v)))
}

source(file.path(root, "scripts", "plot_robustness.R"))   # Figure 5 and Figure S1
cat(sprintf("seed=%d, n_mc=%d per cell, all three tests use the same draws\n", SEED, N_MC))
