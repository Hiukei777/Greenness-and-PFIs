#!/usr/bin/env Rscript
# =============================================================================
# run_all.R -- reproduce every table and figure in the manuscript
#
# Usage (from the repository root):
#   Rscript run_all.R
#
# The steps below are also runnable one at a time, e.g.
#   Rscript R/05_figure2_rcs.R
# because each script sources the prerequisite steps it needs.
# =============================================================================

scripts <- c(
  "R/00_setup.R",
  "R/01_table_s4_models.R",
  "R/02_psm_table1.R",
  "R/03_figure1_balance.R",
  "R/04_table3_conditional_logistic.R",
  "R/05_figure2_rcs.R",
  "R/06_figure3_subgroup.R",
  "R/07_table_s12_quartiles.R",
  "R/08_sensitivity.R",
  "R/09_evalue.R"
)

for (s in scripts) {
  message("\n=== ", s, " ", strrep("=", max(0, 58 - nchar(s))))
  source(s, echo = FALSE)
}

message("\nAll done. Figures and tables were written to the repository root.")
