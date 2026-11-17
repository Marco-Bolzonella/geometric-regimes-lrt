# Figure 7 (Section 5.3): confidence sets for the fifteen coefficients on the
# measured data, from tables/measured_intervals_{1,2}.csv. One row per
# coefficient, ordered by xi(0); for each, the exact, chi2_1 and chibar2
# sets, slightly offset, with the MLE as a dot. Sets of two components are
# drawn as two segments.

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
root <- dirname(this_dir)
source(file.path(root, "core", "plot_theme.R"))

parse_set <- function(s) {
  parts <- regmatches(s, gregexpr("\\[[^]]*\\]", s))[[1]]
  do.call(rbind, lapply(parts, function(p) {
    v <- as.numeric(strsplit(gsub("\\[|\\]", "", p), ",")[[1]])
    data.frame(lo = v[1], hi = v[2])
  }))
}

lev_lab <- c("1" = "68.3% (1σ)", "2" = "95.4% (2σ)")
meth <- c(set_exact = "exact", set_chi2 = "chi2", set_chibar = "chibar")
offs <- c(exact = 0.22, chi2 = 0, chibar = -0.22)

segs <- list(); pts <- list()
for (k in 1:2) {
  d <- read.csv(file.path(root, "tables", sprintf("measured_intervals_%d.csv", k)))
  d <- d[order(d$xi0), ]
  ylab <- sprintf("%s  (ξ = %s)", d$operator, trimws(formatC(signif(d$xi0, 2), format = "fg", digits = 2, big.mark = ",")))
  for (i in seq_len(nrow(d))) for (m in names(meth)) {
    s <- parse_set(d[[m]][i])
    segs[[length(segs) + 1]] <- data.frame(level = lev_lab[[as.character(k)]], row = i, lab = ylab[i],
                                           method = meth[[m]], s)
  }
  pts[[k]] <- data.frame(level = lev_lab[[as.character(k)]], row = seq_len(nrow(d)), chat = d$chat)
  labs_y <- ylab
}
segs <- do.call(rbind, segs); pts <- do.call(rbind, pts)
segs$method <- factor(segs$method, levels = unname(meth))
segs$y <- segs$row + offs[as.character(segs$method)]

p <- ggplot(segs) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = method), linewidth = 0.7) +
  geom_point(data = pts, aes(x = chat, y = row + offs[["exact"]]), size = 0.9) +
  facet_wrap(~level, nrow = 1) +
  scale_y_continuous(breaks = seq_along(labs_y), labels = labs_y, expand = c(0.02, 0.02)) +
  scale_colour_manual(values = c(exact = COL_EXACT, chi2 = COL_CHI2, chibar = COL_CHIBAR),
                      labels = c(exact = "Exact", chi2 = expression(chi[1]^2), chibar = expression(bar(chi)^2)),
                      name = NULL) +
  labs(x = "Wilson coefficient c", y = NULL) +
  theme_paper() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
save_fig(p, file.path(root, "figures", "intervals_forest.png"), height = 4.6)
