# Reproduces the numbers, tables and figures of the paper and of the
# Supplementary Material in order. Run from this folder: Rscript run_all.R
# (a few hours in total).

this_dir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)), fixed = TRUE))
if (this_dir == "") this_dir <- "."

steps <- list(
  c("scripts/validate.R"),                      # self-tests of the geometry and quadrature
  c("scripts/make_geom_schematic.R"),           # Figure 2
  c("scripts/make_regime_figures.R"),           # Figure 3
  c("scripts/make_table1.R"),                   # Table 1, xi and phi
  c("scripts/quadrature_selfcheck.R"),          # Table S1, standard vs reference grid
  c("scripts/quadrature_finesearch.R", "1"),    # Table S1, caption
  c("scripts/quadrature_finesearch.R", "2"),
  c("scripts/tubular_check.R"),                 # Table S1, reference grid vs Eq. (15); Section 5.3 and Appendix F, cQd1 quantiles
  c("scripts/make_coverage_heatmaps.R"),        # Figure 4 and the numbers of Section 4.2
  c("scripts/branch_switching.R"),              # Section 4.2 and Appendix F, branch switching
  c("scripts/make_poisson_robustness.R"),       # Figure 6, Table S2
  c("scripts/make_fig4.R"),                     # Figure S2
  c("scripts/data_summary.R"),                  # Section 5.1 and Appendix F, data and Monte Carlo cost
  c("scripts/measured_intervals.R", "1", "2", "3"), # Section 5.3, Table 1 (MLE, W(0)), Appendix F (3 sigma)
  c("scripts/make_forest_plot.R"),              # Figure 7 (needs measured_intervals_{1,2}.csv)
  c("scripts/make_trajectory_plot.R"),          # Figure 8 (needs measured_intervals_1.csv)
  c("scripts/endpoint_accuracy.R"),             # Table S1, caption, and Section 5.3: interval endpoints
  c("scripts/halfline_check.R"),                # Sections 3.3, 5.3 and Appendix B.4, the half-line
  c("scripts/partial_correction.R"),            # Appendix F, partial corrections
  c("scripts/coverage_length.R"),               # Section 5.4 and Appendix F, coverage and length
  c("scripts/section4_numbers.R"),              # remaining numbers of Sections 4-5 and Appendix F
  c("scripts/bernlochner_check.R")              # Appendix D
)

run_step <- function(step) {
  cat(sprintf("\n%s\n%s\n%s\n", strrep("=", 70), paste(step, collapse = " "), strrep("=", 70)))
  status <- system2("Rscript", c(file.path(this_dir, step[1]), step[-1]))
  if (status != 0) {
    cat(sprintf("FAILED at %s -- stopping.\n", step[1]))
    quit(status = 1)
  }
}
for (step in steps) run_step(step)

# Figures 5 and S1: the Monte Carlo grid is computed in chunks; rerun until COMPLETE
cat(sprintf("\n%s\nscripts/make_robustness_heatmaps.R\n%s\n", strrep("=", 70), strrep("=", 70)))
repeat {
  out <- system2("Rscript", file.path(this_dir, "scripts/make_robustness_heatmaps.R"),
                 stdout = TRUE, stderr = TRUE)
  cat(paste(out, collapse = "\n"), "\n")
  if (any(grepl("^FAILED|Error", out))) {
    cat("FAILED at scripts/make_robustness_heatmaps.R -- stopping.\n")
    quit(status = 1)
  }
  if (any(grepl("^COMPLETE", out))) break
}

# every number quoted in the text that can be recomputed from the saved results
run_step("scripts/check_draft_numbers.R")
cat("\nAll steps completed. See figures/ and tables/.\n")
