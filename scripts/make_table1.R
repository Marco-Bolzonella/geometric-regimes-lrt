# Table 1 (Section 5.2): (xi, phi, regime) of the 15 coefficients.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "load_operators.R"))

df <- build_table1(file.path(root, "data"))
out <- file.path(root, "tables", "table1_xi_phi.csv")
df$xi <- round(df$xi, 3)
df$phi_deg <- round(df$phi_deg, 2)
write.csv(df[, c("operator", "xi", "phi_deg", "regime")], out, row.names = FALSE)
cat(sprintf("wrote %s (%d operators)\n", out, nrow(df)))
