# =============================================================================
# Protocol_Install.R
#
# Installs every R package used by the protocol (Protocol_Combined.qmd) and its
# scripts, then prints the versions installed.
#
# Before running this script:
#   1. Install R and RStudio. Recent versions of RStudio include Quarto, which
#      is needed to render the protocol.
#   2. Install JAGS (version 4.3.1 or newer) from https://mcmc-jags.sourceforge.io/
#      JAGS is a separate program, not an R package. RoBMA needs it.
#
# Run time: a few minutes, depending on how many packages are already installed.
# =============================================================================

# ---- Step 1: packages from CRAN --------------------------------------------
# Only the packages that are not yet installed are installed.
cran_packages <- c(
  "here",          # resolves file paths relative to the project root (.Rproj)
  "metafor",       # meta-analytic models
  "clubSandwich",  # cluster-robust (CR2) standard errors
  "tidyverse",     # data manipulation and plotting
  "haven",         # reading Stata .dta files
  "readxl",        # reading Excel output in the display chunks
  "writexl",       # writing Excel output
  "openxlsx",      # writing Excel output with formatting (study weights table)
  "knitr",         # tables and figures in the display chunks
  "future",        # parallel processing (selection-model bootstrap)
  "progressr",     # bootstrap progress display
  "metaselection", # step-function selection models
  "posterior"      # posterior draws for RoBMA predictions
)

missing <- cran_packages[!cran_packages %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

# ---- Step 2: RoBMA and BayesTools ------------------------------------------
# The protocol fits RoBMA with the CHE model using RoBMA.mv() (Section 5.4),
# which first appeared in RoBMA 4.1.1. Both packages come from CRAN. They are
# handled separately from Step 1 because an older RoBMA (for example 4.0.0)
# may already be installed: Step 1 would skip it, but it must be updated.
# This step is skipped if RoBMA 4.1.1 or later is already installed.
if (!requireNamespace("RoBMA", quietly = TRUE) ||
    packageVersion("RoBMA") < "4.1.1") {
  install.packages(c("BayesTools", "RoBMA"))
}

# ---- Step 3: report the versions installed ---------------------------------
# The results reported in the protocol were produced with RoBMA 4.1.1,
# BayesTools 0.4.0, and metaselection 0.1.5.
cat("\nInstalled versions:\n")
for (p in c("RoBMA", "BayesTools", "metaselection", "metafor", "clubSandwich")) {
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))
}
cat("\nIf RoBMA was installed in this session, restart R before running the scripts.\n")
