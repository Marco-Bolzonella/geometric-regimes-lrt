# Supplementary Material, Appendix D: checks that the (Delta0, Delta1) of
# Bernlochner et al. (2023), written in c, equal ours times ||T||^6 and
# ||T||^9, and that their standard deviations are 12 Nbar^3 and
# 108 (1 - rho^2)^(1/2) Lbar Nbar^4.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "load_operators.R"))
set.seed(20260909)
raw <- load_raw(file.path(root, "data")); M <- whitening_matrix(raw$cov)
ops <- build_table1(file.path(root, "data"))$operator
for (op in ops) for (c0 in c(0, 0.3)) {
  l <- as.numeric(M %*% sapply(raw$bins, function(b) raw$param[[b]][[op]]))
  n <- as.numeric(M %*% sapply(raw$bins, function(b) raw$param[[b]][[paste0(op, "*", op)]]))
  L2 <- sum(l^2); M2 <- sum(l * n); N2 <- sum(n^2); rho <- M2 / sqrt(L2 * N2)
  T_ <- l + 2 * c0 * n; Nv <- 2 * n; nT <- sqrt(sum(T_^2))
  xi <- sqrt(sum(Nv^2)) / nT^2; e1 <- T_ / nT; a <- sum(Nv * e1); e2 <- (Nv - a * e1); e2 <- e2 / sqrt(sum(e2^2))
  ph <- atan2(sum(Nv * e2), a)
  res <- replicate(4000, {
    z <- rnorm(length(l)); x <- c0 * l + c0^2 * n + z
    Xl <- sum(x * l); Xn <- sum(x * n)
    # their cubic: 2 N2 c^3 + 3 M2 c^2 + (L2 - 2 Xn) c - Xl = 0
    Ab <- 2 * N2; Bb <- 3 * M2; Cb <- L2 - 2 * Xn; Db <- -Xl
    d0b <- Bb^2 - 3 * Ab * Cb; d1b <- 2 * Bb^3 - 9 * Ab * Bb * Cb + 27 * Ab^2 * Db
    zu <- sum(z * e1); zv <- sum(z * e2)
    A <- 0.5 * xi^2; B <- 1.5 * xi * cos(ph); C <- 1 - xi * cos(ph) * zu - xi * sin(ph) * zv; D <- -zu
    d0 <- B^2 - 3 * A * C; d1 <- 2 * B^3 - 9 * A * B * C + 27 * A^2 * D
    c(d0b / d0, d1b / d1, d0b, d1b)
  })
  cat(sprintf("%-6s c0=%.1f: ratio Delta0 = %.6f (|T|^6 = %.6f, sd %.2e); ratio Delta1 = %.6f (|T|^9 = %.6f, sd %.2e)\n",
              op, c0, mean(res[1, ]), nT^6, sd(res[1, ]), mean(res[2, ]), nT^9, sd(res[2, ])))
  cat(sprintf("         sd(Delta0) = %.4g vs 12 N^3 = %.4g; sd(Delta1) = %.4g vs 108 sqrt(1-rho^2) L N^4 = %.4g; cor = %.3f\n",
              sd(res[3, ]), 12 * N2^1.5, sd(res[4, ]), 108 * sqrt(1 - rho^2) * sqrt(L2) * N2^2, cor(res[3, ], res[4, ])))
}
