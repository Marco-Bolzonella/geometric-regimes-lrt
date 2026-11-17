# Code and data for "Beyond the Standard Model of Particle Physics: Geometric Regimes of the Likelihood Ratio Statistic"

This package reproduces every number, table and figure of the paper and of
its Supplementary Material, except Figure 1, which is a schematic
illustration.

## Requirements

R with the packages `yaml` and `ggplot2` (version 3.4 or later; `scales` is
installed with it). The results were produced with R 4.6.0, ggplot2 4.0.3
and yaml 2.3.12. Figures use Unicode Greek letters, so R must run in a UTF-8
locale (`core/plot_theme.R` switches to `C.UTF-8` if needed).

## Running

From this folder:

    Rscript run_all.R

runs the self-tests and then every script in order. Scripts must be run
with `Rscript`, which they use to locate the package folder. The full run takes a
few hours on a laptop, mostly for the confidence intervals, which use the
reference quadrature grid, and for the Monte Carlo studies. Any script can
be run on its own, from any folder, e.g. `Rscript scripts/make_fig4.R`; outputs go to
`tables/` and `figures/`, which already contain the results reported in the
paper.

## Which script makes what

Numbers refer to the paper; tables and figures numbered S1, S2, ... are in
the Supplementary Material.

| Item | Script | Output |
|---|---|---|
| Figure 2 (Section 2.3), geometry of the model curve | `scripts/make_geom_schematic.R` | `figures/geom_schematic.pdf` |
| Figure 3 (Section 3), the three regimes | `scripts/make_regime_figures.R` | `figures/regimes.png` |
| Section 4.1 and Table S1, accuracy of the quadrature; Section 5.3 and Appendix F, exact quantiles at the endpoints and for cQd1 | `scripts/quadrature_selfcheck.R`, `scripts/quadrature_finesearch.R` (argument `1`, then `2`), `scripts/tubular_check.R`, `scripts/endpoint_accuracy.R` | `tables/quadrature_selfcheck.csv`, `tables/quadrature_finesearch_{1,2}.csv`, `tables/tubular_check.csv`, `tables/endpoint_accuracy.csv` |
| Figure 4, Section 4.2 and Appendix F (higher levels), error rates of the chi2_1 and chibar2 tests | `scripts/make_coverage_heatmaps.R` | `figures/error_rate_heatmaps.png`, `tables/coverage_grids_levels.rds` |
| Section 4.2 and Appendix F, branch switching | `scripts/branch_switching.R` | `tables/branch_switching.csv` |
| Figures 5 and S1, Section 4.3, heavy tails | `scripts/make_robustness_heatmaps.R` (computes; rerun until it prints `COMPLETE`), `scripts/plot_robustness.R` (draws) | `tables/robustness_grid_levels.rds`, `figures/robustness_heatmaps.png`, `figures/robustness_heatmaps_all.png` |
| Figure 6, Table S2, Section 4.4 (including the equivalent counts), Poisson counts | `scripts/make_poisson_robustness.R` (computes), `scripts/plot_poisson.R` (draws) | `tables/poisson_robustness.csv`, `figures/poisson_size.png` |
| Section 5.1 and Appendix F, measured data and cost of Monte Carlo calibration | `scripts/data_summary.R` | console |
| Sections 3.3 and 5.3 and Appendix B.4, the half-line case sin(phi) = 0 | `scripts/halfline_check.R` | `tables/halfline_check.csv` |
| Table 1, xi and phi at c0 = 0 | `scripts/make_table1.R` | `tables/table1_xi_phi.csv` |
| Table 1 (MLE, W(0)), Section 5.3, Appendix F (3 sigma) | `scripts/measured_intervals.R` (arguments `1 2 3`) | `tables/measured_intervals_{1,2,3}.csv` |
| Figure 7, confidence sets | `scripts/make_forest_plot.R` | `figures/intervals_forest.png` |
| Figure 8 and Section 5.3, geometry along the confidence sets | `scripts/make_trajectory_plot.R` | `figures/trajectory_plot.png`, `tables/partial_correction.csv` |
| Section 5.4 and Appendix F, coverage, length and number of components of the confidence sets | `scripts/coverage_length.R` | `tables/coverage_length.csv` |
| Appendix F, partial corrections | `scripts/partial_correction.R` | `tables/partial_correction_all.csv` |
| Figure S2, exact and simulated quantiles | `scripts/make_fig4.R` | `figures/fig4_qqplot.png` |
| Remaining numbers of Sections 4 and 5 and Appendix F | `scripts/section4_numbers.R` | console |
| Appendix D, equivalence with Bernlochner et al. (2023) | `scripts/bernlochner_check.R` | console |
| Numbers in the text recomputed from the saved results | `scripts/check_draft_numbers.R` | console |

`scripts/validate.R` runs the self-tests of the geometry and of the
quadrature.

## Core code

- `core/geom_core.R`: the likelihood ratio statistic W(c0) for given
  (z_u, z_v, xi, phi), from the real roots of the cubic stationarity
  equation, choosing the global minimiser between the two outer roots.
- `core/geom_quantiles.R`: tail probabilities and quantiles of W(c0) by
  two-dimensional quadrature over the Gaussian signal plane, on the
  standard grid (48 angular x 400 radial nodes) or on the reference grid
  (192 x 800).
- `core/geom_tubular.R`: the exact density of W(c0) in tubular coordinates,
  Eq. (15), by adaptive quadrature in tau; used to check the plane
  quadrature.
- `core/load_operators.R`: reads the inputs, whitens them and computes
  (xi, phi).
- `core/geom_poisson.R`: Poisson simulation, with the likelihood ratio
  computed from the Gaussian working likelihood.
- `core/plot_theme.R`: figure size, fonts and labels.

## Data

The inputs of the case study come from the HEPData record of the ATLAS
measurement (JHEP 06 (2022) 063), https://www.hepdata.net/record/ins2037744
(doi:10.17182/hepdata.134011.v2).

- `data/hepdata_134011_table2_abs_ptthad.yaml` (included): Table 2 of the
  record, the absolute differential cross-section in the transverse momentum
  of the hadronically decaying top quark and the Powheg+Pythia8 prediction.
  The data used in the paper are the ratios of the two.
- `data/covariance_ttbar_boosted_exp_v1_0.yaml` (not included): the 8 x 8
  covariance matrix of the ratios, as a YAML list of rows under the key
  `covariance`.
- `data/param_ttbar_boosted_quad_v1_0.yaml` (not included): the bin-wise
  linear and quadratic sensitivities (l_i, q_i) of the fifteen Wilson
  coefficients, under the key `parameterisation`, with one entry per bin
  and keys `<op>` (l_i) and `<op>*<op>` (q_i).

The last two files are internal inputs of the ATLAS analysis and are not
distributed; the information they contain can be derived from the HEPData
record above. Scripts that need them stop with a message if they are
absent. The tables and figures in `tables/` and `figures/` are those
obtained with them.

## Reproducibility

Every script that draws random numbers is seeded, with seed 20260909 (offset
per row in `make_robustness_heatmaps.R`) except for the illustrative points
of Figure 3, which use seed 123, so reruns give identical results.
`make_coverage_heatmaps.R` reuses `tables/coverage_grids_levels.rds` if it
exists, and so only redraws Figure 4; delete that file to recompute the grid. The
quadrature results are deterministic. The figures are 6.5 inches wide at
300 dpi, except Figure 2, which is a vector PDF.

## Licence

The code is released under the MIT licence (`LICENSE`). The HEPData record
(https://www.hepdata.net/record/ins2037744) is distributed under its own
terms (CC0).
