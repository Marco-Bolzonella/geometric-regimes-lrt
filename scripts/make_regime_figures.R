# Figure 3 (Section 3): the curve for xi = 0.1, 6 and 500 (phi = pi/2), 50 draws
# of Z_L and their projections onto the curve, one panel per regime.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "plot_theme.R"))

XI <- c(0.1, 6, 500)
SEED <- 123
POINTS <- 50

panel_data <- function(xi) {
  lab <- sprintf("ξ = %s", format(xi))
  u <- seq(-4, 4, length.out = 4000)
  curve <- data.frame(panel = lab, u = u, v = 0.5 * xi * u^2)
  set.seed(SEED)
  pts <- data.frame(panel = lab, zu = rnorm(POINTS), zv = rnorm(POINTS))
  # nearest point on the curve
  k <- vapply(seq_len(POINTS), function(i)
    which.min((curve$u - pts$zu[i])^2 + (curve$v - pts$zv[i])^2), integer(1))
  pts$u_hat <- curve$u[k]; pts$v_hat <- curve$v[k]
  list(curve = curve, pts = pts)
}
d <- lapply(XI, panel_data)
curve <- do.call(rbind, lapply(d, `[[`, "curve"))
pts <- do.call(rbind, lapply(d, `[[`, "pts"))
lev <- unique(curve$panel)
curve$panel <- factor(curve$panel, levels = lev); pts$panel <- factor(pts$panel, levels = lev)

p <- ggplot() +
  geom_segment(data = pts, aes(zu, zv, xend = u_hat, yend = v_hat), colour = "grey70", linewidth = 0.3) +
  geom_line(data = curve, aes(u, v), colour = COL_CURVE, linewidth = 0.8) +
  geom_point(data = pts, aes(zu, zv), colour = "#D55E00", size = 0.9) +
  annotate("point", x = 0, y = 0, size = 1.8) +
  facet_wrap(~panel, nrow = 1) +
  coord_fixed(xlim = c(-2, 2), ylim = c(-2, 2.5), expand = FALSE) +
  scale_x_continuous(breaks = -1:1) +
  labs(x = expression(z[u]), y = expression(z[v])) +
  theme_paper()
save_fig(p, file.path(root, "figures", "regimes.png"), height = 2.6)
