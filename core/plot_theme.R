# Figure settings: 6.5 in wide at 300 dpi, Unicode Greek letters (so that
# phi prints as the varphi of the paper).

library(ggplot2)

if (!l10n_info()$`UTF-8`) invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

FIG_WIDTH <- 6.5      # inches; scaled down to \textwidth in LaTeX
FIG_DPI   <- 300
BASE_SIZE <- 9

theme_paper <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(strip.background = element_rect(fill = "grey92", colour = "grey40"),
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.9, "lines"),
          plot.margin = margin(4, 14, 4, 4))
}

# one colour per object across all figures (Okabe-Ito); the regimes get the
# three Okabe-Ito colours not used for the quantiles
COL_EXACT   <- "black"
COL_CHI2    <- "#D55E00"
COL_CHIBAR  <- "#009E73"
COL_CURVE   <- "#0072B2"
COL_REGIME  <- c("Regime I" = "#56B4E9", "Regime II" = "#E69F00", "Regime III" = "#CC79A7")

LAB_XI      <- "ξ (log scale)"
LAB_PHI     <- "φ (degrees)"
DEG         <- "°"

# shared pieces of the heatmaps over the (xi, phi) plane (Figures 4, 5, S1)
scale_x_xi <- function()
  scale_x_log10(expand = c(0, 0), breaks = 10^(-1:4), labels = parse(text = paste0("10^", -1:4)))
scale_y_phi <- function() scale_y_continuous(expand = c(0, 0), breaks = seq(0, 180, 45))
geom_operators <- function(ops)
  geom_point(data = ops, aes(xi, phi_deg), inherit.aes = FALSE,
             shape = 21, fill = "gold", colour = "black", size = 0.9, stroke = 0.3)
theme_heatmap <- function()
  theme_paper() + theme(legend.key.width = unit(0.6, "lines"), legend.key.height = unit(1.4, "lines"),
                        panel.spacing.x = unit(0.8, "lines"))

save_fig <- function(p, path, height, width = FIG_WIDTH) {
  ggsave(path, p, width = width, height = height, dpi = FIG_DPI, bg = "white")
  cat(sprintf("wrote %s (%.1f x %.1f in, %d dpi)\n", path, width, height, FIG_DPI))
}
