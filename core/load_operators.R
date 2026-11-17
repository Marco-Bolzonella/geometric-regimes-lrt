# Reads the sensitivities and covariance of the case study (8 bins, 15
# Wilson coefficients) and computes (xi, phi) for each coefficient.

library(yaml)

load_raw <- function(data_dir_path) {
  needed <- c("param_ttbar_boosted_quad_v1_0.yaml", "covariance_ttbar_boosted_exp_v1_0.yaml")
  missing <- needed[!file.exists(file.path(data_dir_path, needed))]
  if (length(missing))
    stop("missing input file(s) in data/: ", paste(missing, collapse = ", "),
         "\nThese are not distributed; see the Data section of README.md.", call. = FALSE)
  param <- yaml.load_file(file.path(data_dir_path, "param_ttbar_boosted_quad_v1_0.yaml"))$parameterisation
  cov_list <- yaml.load_file(file.path(data_dir_path, "covariance_ttbar_boosted_exp_v1_0.yaml"))$covariance
  cov <- do.call(rbind, lapply(cov_list, unlist))
  bins <- names(param)
  stopifnot(nrow(cov) == length(bins), ncol(cov) == length(bins))
  list(param = param, bins = bins, cov = cov)
}

# whitening matrix M with t(M) %*% M = solve(Sigma0)
whitening_matrix <- function(cov) {
  D <- solve(cov)
  Lc <- chol(D)
  Lc
}

# (xi, phi) at c0 = 0: T = M l, N = 2 M q
xi_phi_for_operator <- function(op, param, bins, M) {
  l <- sapply(bins, function(b) param[[b]][[op]])
  q <- sapply(bins, function(b) param[[b]][[paste0(op, "*", op)]])
  T_ <- as.numeric(M %*% l)
  N_ <- as.numeric(2 * (M %*% q))
  normT <- sqrt(sum(T_^2))
  e1 <- T_ / normT
  a <- sum(N_ * e1)
  Nperp <- N_ - a * e1
  normN <- sqrt(sum(N_^2))
  xi <- normN / normT^2
  phi <- atan2(sqrt(sum(Nperp^2)), a)  # radians, in [0, pi]
  c(xi = xi, phi = phi)
}

ALL_OPERATORS <- c(
  "ctGRe", "ctj8", "cQj18", "cQu8", "ctu8", "cQd8", "ctd8",
  "cQj38", "cQj31", "ctu1", "cQj11", "cG", "ctj1", "cQu1", "cQd1"
)

regime <- function(xi) {
  ifelse(xi <= 5, "I", ifelse(xi <= 200, "II", "III"))
}

build_table1 <- function(data_dir_path) {
  raw <- load_raw(data_dir_path)
  M <- whitening_matrix(raw$cov)
  rows <- lapply(ALL_OPERATORS, function(op) {
    xp <- xi_phi_for_operator(op, raw$param, raw$bins, M)
    data.frame(operator = op, xi = xp["xi"], phi_rad = xp["phi"],
               phi_deg = xp["phi"] * 180 / pi,
               regime = regime(xp["xi"]), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df[order(df$xi), ]
}

if (sys.nframe() == 0) {
  # Rscript core/load_operators.R prints Table 1
  this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
  ddir <- file.path(dirname(this_dir), "data")
  df <- build_table1(ddir)
  cat(sprintf("n_operators = %d\n", nrow(df)))
  print(df, row.names = FALSE)
  cat("\nRegime counts:\n")
  print(table(df$regime))
}
