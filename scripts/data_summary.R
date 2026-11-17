# Section 5.1, description of the measured data, and Appendix F, the Monte Carlo
# cost arithmetic of the same section. Prints to the console.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
library(yaml)

hep <- yaml.load_file(file.path(root, "data", "hepdata_134011_table2_abs_ptthad.yaml"))
meas <- sapply(hep$dependent_variables[[1]]$values, function(v) v$value)
pred <- sapply(hep$dependent_variables[[2]]$values, function(v) v$value)
cov_file <- file.path(root, "data", "covariance_ttbar_boosted_exp_v1_0.yaml")
if (!file.exists(cov_file))
  stop("missing input file data/covariance_ttbar_boosted_exp_v1_0.yaml\n",
       "It is not distributed; see the Data section of README.md.", call. = FALSE)
S <- do.call(rbind, lapply(yaml.load_file(cov_file)$covariance, unlist))
D <- solve(S)

Y <- meas / pred
r <- Y - 1
cat(sprintf("ratios Y_i: %.3f to %.3f\n", min(Y), max(Y)))
cat(sprintf("chi2 of the Standard Model: %.2f on %d degrees of freedom\n", drop(t(r) %*% D %*% r), length(Y)))
cat(sprintf("standard deviations of the ratios: %.1f%% to %.1f%%\n", 100 * min(sqrt(diag(S))), 100 * max(sqrt(diag(S)))))
R <- cov2cor(S)
cat(sprintf("largest correlation: %.3f\n", max(R[upper.tri(R)])))
e <- eigen(S, symmetric = TRUE)
cat(sprintf("leading principal component: %.1f%% of the variance, loadings of one sign: %s\n",
            100 * e$values[1] / sum(e$values), all(sign(e$vectors[, 1]) == sign(e$vectors[1, 1]))))
one <- rep(1, length(Y))
k <- drop(t(one) %*% D %*% Y / (t(one) %*% D %*% one))
cat(sprintf("common rescaling of the prediction: %.3f, leaving chi2 = %.2f\n", k, drop(t(Y - k) %*% D %*% (Y - k))))

# Monte Carlo cost: standard error of an empirical p-quantile, {p(1-p)/n}^{1/2}/f(w_p),
# at the chi2_1 density f(1); accuracies of the standard (0.167 points) and
# reference (0.014 points) grids at 1 sigma (Table S1); 300 scanned null values per interval
p <- 2 * pnorm(1) - 1
f1 <- dchisq(1, 1)
for (acc in c(0.167, 0.014)) {
  dw <- acc / 100 / f1
  cat(sprintf("accuracy %.3f points: %.4f in the threshold, about %.1e draws\n", acc, dw, p * (1 - p) / (acc / 100)^2))
}
cat(sprintf("fits for one interval at the standard-grid accuracy: %.1e\n", round(p * (1 - p) / (0.167 / 100)^2, -3) * 300))
