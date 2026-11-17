# Recomputes from the saved results every number quoted in the paper and the
# Supplementary Material that can be checked without new simulation; prints one line per claim.
this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
setwd(dirname(this_dir))
T <- function(...) cat(sprintf(...), "\n")
tb <- function(f) read.csv(file.path("tables", f))

T("--- 4.1 quadrature")
T("finesearch refined: 1s %.4f, 2s %.4f (prob)", tb("quadrature_finesearch_1.csv")$refined_max_pp/100, tb("quadrature_finesearch_2.csv")$refined_max_pp/100)
tc <- tb("tubular_check.csv"); ea <- tb("endpoint_accuracy.csv")
T("ref vs tubular max %.5f; endpoints max %.5f", max(abs(tc$diff_ref_pp))/100, max(abs(ea$diff_pp))/100)
qs <- tb("quadrature_selfcheck.csv"); print(qs[, c("sigma", "max_abs_error_pp")])
print(aggregate(abs(diff_ref_pp) ~ k_sigma, tc, max))

T("--- 4.2 error-rate grid")
G <- readRDS("tables/coverage_grids_levels.rds")
for (k in 1:3) for (m in c("chi2", "chernoff")) {
  g <- G[[k]][[m]]; r <- 1 - g$achieved; a <- 1 - G[[k]]$cov_nom
  T("%ds %-8s rate [%.3f, %.3f] ratio max %.2f min %.2f", k, m, min(r), max(r), max(r/a), min(r/a))
}
g <- G[[1]]$chernoff; xs <- unique(g$xi); x3 <- xs[which.min(abs(log(xs) - log(1000)))]
T("chibar 1s at xi=%.0f: rate [%.3f, %.3f]", x3, min(1 - g$achieved[g$xi == x3]), max(1 - g$achieved[g$xi == x3]))
g <- G[[1]]$chi2; a <- 1 - G[[1]]$cov_nom
cross <- sapply(split(g, g$phi_deg), function(d) { d <- d[order(d$xi), ]; d$xi[which(1 - d$achieved < a)[1]] })
T("chi2 1s: xi where rate first drops below nominal, range over phi [%.2f, %.2f]", min(cross), max(cross))
p90 <- unique(g$phi_deg)[which.min(abs(unique(g$phi_deg) - 90))]
s <- g$phi_deg == p90 & g$xi >= 0.5 & g$xi <= 1.5
T("chi2 1s at phi=%.1f, xi in [0.5,1.5]: rate [%.4f, %.4f]", p90, min(1 - g$achieved[s]), max(1 - g$achieved[s]))
bs <- tb("branch_switching.csv"); print(bs)

T("--- 4.3 robustness")
R <- readRDS("tables/robustness_grid_levels.rds")
for (k in 1:2) for (s in names(R$results[[k]])) for (m in c("exact", "chi2", "chibar")) {
  v <- R$results[[k]][[s]][[m]]
  T("%ds %-6s %-6s [%.4f, %.4f]  large-xi(>1000) min %.4f", k, s, m, min(v), max(v), min(v[, R$xi_grid > 1000]))
}
for (k in 1:2) { a <- R$alphas[k]
  wn <- max(sapply(c("t5", "t3"), function(s) max(abs(R$results[[k]][[s]]$exact - a))))
  wc <- max(abs(R$results[[k]]$Normal$chi2 - a))
  T("%ds worst |exact-nom| heavy tails %.4f vs chi2 normal %.4f", k, wn, wc) }

T("--- 5.3 intervals")
for (k in 1:3) { d <- tb(sprintf("measured_intervals_%d.csv", k))
  reg1 <- d$xi0 <= 5
  T("%ds pct chibar [%.1f, %.1f]; chi2 regI [%.1f, %.1f]; chi2 other [%.1f, %.1f]", k, min(d$pct_chibar), max(d$pct_chibar),
    min(d$pct_chi2[reg1]), max(d$pct_chi2[reg1]), min(d$pct_chi2[!reg1]), max(d$pct_chi2[!reg1]))
  has0 <- function(s) sapply(s, function(x) { p <- regmatches(x, gregexpr("-?[0-9.]+", x))[[1]]; v <- as.numeric(p)
    any(v[c(TRUE, FALSE)] <= 0 & v[c(FALSE, TRUE)] >= 0) })
  T("%ds sets without 0: exact {%s} chi2 {%s} chibar {%s}", k, paste(d$operator[!has0(d$set_exact)], collapse = ","),
    paste(d$operator[!has0(d$set_chi2)], collapse = ","), paste(d$operator[!has0(d$set_chibar)], collapse = ","))
}
d <- tb("measured_intervals_1.csv"); print(d[d$operator %in% c("ctj8", "cQu8"), c("operator", "set_exact", "set_chi2")])
T("Table 1 xi gaps: sorted xi %s", paste(signif(sort(d$xi0), 3), collapse = " "))

T("--- 5.4 coverage")
cv <- tb("coverage_length.csv"); se <- sqrt(cv$nominal * (1 - cv$nominal) / cv$B)
T("exact max |dev|/se %.2f; se 1s %.4f 2s %.4f", max(abs(cv$cover_exact - cv$nominal) / se), se[cv$level_sigma == 1][1], se[cv$level_sigma == 2][1])
for (k in 1:2) { s <- cv$level_sigma == k
  T("%ds chi2 [%.4f, %.4f] chibar [%.4f, %.4f]", k, min(cv$cover_chi2[s]), max(cv$cover_chi2[s]), min(cv$cover_chibar[s]), max(cv$cover_chibar[s])) }
x0 <- setNames(d$xi0, d$operator); cv$xi0 <- x0[cv$operator]
s <- cv$level_sigma == 1
T("1s chi2 under-covers for: %s", paste(unique(cv$operator[s & cv$cover_chi2 < cv$nominal - 3 * se]), collapse = ","))
T("1s chi2 over-covers for: %s", paste(unique(cv$operator[s & cv$cover_chi2 > cv$nominal + 3 * se]), collapse = ","))
for (k in 1:2) { s <- cv$level_sigma == k
  T("%ds multi-component exact: regII/III [%.3f, %.3f]  regI max %.3f; fixed regII/III chi2 [%.3f,%.3f] chibar [%.3f,%.3f]", k,
    min(cv$frac_multi_component[s & cv$xi0 > 5]), max(cv$frac_multi_component[s & cv$xi0 > 5]), max(cv$frac_multi_component[s & cv$xi0 <= 5]),
    min(cv$frac_multi_chi2[s & cv$xi0 > 5]), max(cv$frac_multi_chi2[s & cv$xi0 > 5]),
    min(cv$frac_multi_chibar[s & cv$xi0 > 5]), max(cv$frac_multi_chibar[s & cv$xi0 > 5])) }
cv$r1 <- cv$mean_len_chi2 / cv$mean_len_exact; cv$rb <- cv$mean_len_chibar / cv$mean_len_exact
for (k in 1:2) { s <- cv$level_sigma == k
  T("%ds len ratio chibar [%.2f, %.2f]; chi2 regI [%.2f, %.2f] regII/III [%.2f, %.2f]", k, min(cv$rb[s]), max(cv$rb[s]),
    min(cv$r1[s & cv$xi0 <= 5]), max(cv$r1[s & cv$xi0 <= 5]), min(cv$r1[s & cv$xi0 > 5]), max(cv$r1[s & cv$xi0 > 5])) }
T("max |shift_chi2 - shift_exact| / len %.3f; |shift_chibar - shift_exact| / len %.3f",
  max(abs(cv$shift_chi2 - cv$shift_exact) / cv$mean_len_exact), max(abs(cv$shift_chibar - cv$shift_exact) / cv$mean_len_exact))
