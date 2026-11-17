# Figure S2 (Appendix F): exact, chi2_1 and chibar2 quantiles of W(c0) against
# 3e5 draws under the working model, one coefficient per regime.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "geom_core.R"))
source(file.path(root, "core", "geom_quantiles.R"))
source(file.path(root, "core", "load_operators.R"))
source(file.path(root, "core", "plot_theme.R"))

SEED    <- 20260909
N_MC    <- 300000
N_PROBS <- 40
REPRESENTATIVE <- c(ctGRe = "Regime I", cQj31 = "Regime II", cQd1 = "Regime III")

operators <- build_table1(file.path(root, "data"))
rownames(operators) <- operators$operator
probs <- seq(0.01, 0.99, length.out = N_PROBS)

set.seed(SEED)
rows <- list()
for (op in names(REPRESENTATIVE)) {
  xi  <- operators[op, "xi"]
  phi <- operators[op, "phi_rad"]
  Wmc   <- compute_w_vectorized(rnorm(N_MC), rnorm(N_MC), xi, phi)
  emp_q <- quantile(Wmc, probs, names = FALSE)
  theo_q   <- vapply(probs, function(p) quantile_W(xi, phi, 1 - p), numeric(1))
  chi2_q   <- qchisq(probs, df = 1)
  chibar_q <- qchisq(pmax(2 * probs - 1, 0), df = 1)  # 0 for probs <= 1/2 (atom at zero)
  lab <- sprintf("%s (%s)\nξ = %s, φ = %.0f%s", op, REPRESENTATIVE[[op]],
                 trimws(formatC(signif(xi, 2), format = "fg", digits = 2, big.mark = ",")), phi * 180 / pi, DEG)
  rows[[op]] <- rbind(
    data.frame(panel = lab, theo = theo_q,   emp = emp_q, method = "exact"),
    data.frame(panel = lab, theo = chi2_q,   emp = emp_q, method = "wilks"),
    data.frame(panel = lab, theo = chibar_q, emp = emp_q, method = "chernoff"))
}
df <- do.call(rbind, rows)
df$panel  <- factor(df$panel, levels = unique(df$panel))
df$method <- factor(df$method, levels = c("exact", "wilks", "chernoff"))

meth_lab <- c(exact = "Exact", wilks = expression(chi[1]^2), chernoff = expression(bar(chi)^2))
p <- ggplot(df, aes(theo, emp, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
  geom_point(size = 1.1) +
  facet_wrap(~panel, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(exact = COL_EXACT, wilks = COL_CHI2, chernoff = COL_CHIBAR),
                      labels = meth_lab, name = NULL) +
  scale_shape_manual(values = c(exact = 16, wilks = 4, chernoff = 17),
                     labels = meth_lab, name = NULL) +
  labs(x = "Theoretical quantile", y = "Monte Carlo quantile") +
  theme_paper() + theme(legend.position = "bottom")

save_fig(p, file.path(root, "figures", "fig4_qqplot.png"), height = 2.9)
cat(sprintf("seed=%d, n_mc=%d per Wilson coefficient\n", SEED, N_MC))
