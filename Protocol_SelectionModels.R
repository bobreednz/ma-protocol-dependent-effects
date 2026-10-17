# =============================================================================
# Protocol_SelectionModels.R
# =============================================================================
# Purpose: intercept-only step-function selection models (Section 5.3).
#
#   3PSM: one threshold at one-sided p = 0.025 (two-sided significance at 5%).
#         lambda1 = relative reporting probability of everything that is not
#         significant and positive (non-significant positive AND wrong-signed).
#   4PSM: thresholds at one-sided p = 0.025 and 0.5. Splits the 3PSM category
#         in two:
#           lambda1 = positive but not significant (0.025 < p <= 0.5)
#           lambda2 = wrong sign (p > 0.5)
#         4PSM contains 3PSM as a special case (lambda2 = lambda1). Comparing
#         lambda1 and lambda2 shows whether selection operates on sign; the
#         protocol (Section 6) uses that comparison, with the z-plot, to decide
#         which of the two models to rely on.
#
# A lambda of 1 means no selection: results in that category are as likely to
# be reported as significant positive results. A lambda well below 1 means
# results in that category are under-reported.
#
# Estimation: penalized marginal likelihood (PML) with the weak default priors
# of Pustejovsky, Citkowicz, and Joshi (2026). Dependence among estimates from
# the same study is handled through a two-stage cluster bootstrap (R = 1999),
# which resamples studies and then estimates within studies.
#
# Note on the estimand: these models treat each estimate as marginally
# independent with a single heterogeneity parameter (tau), as in a standard
# random-effects model. Their mean is therefore an average across estimates,
# whereas the CHE mean gives each study more equal weight.
#
# This script replaces the selection-model part of Protocol_AppendixPubBias.R.
# The models, priors and seeds are identical, so the estimates are unchanged.
#
# Output:
#   Table_SelectionModels.xlsx   Sheet "Display": formatted table
#                                Sheet "Values":  numeric values used by
#                                                 Protocol_Diagnostics.R
#   sel_3psm_intercept.rds, sel_4psm_intercept.rds   fitted models
# Run time: with refit <- FALSE (the default), a few seconds; the script loads
# the saved fits. With refit <- TRUE, about 3 minutes with four cores, longer
# with fewer (two bootstrap fits).
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

# ---- Settings ------------------------------------------------------------------
refit <- FALSE   # TRUE refits both models with the bootstrap (slow); FALSE loads saved fits

# ---- Packages ----------------------------------------------------------------
# metaselection is on CRAN (version 0.3.0 or later). The results reported in
# the protocol were produced with development version 0.1.5.
packages <- c("here", "metaselection", "future", "progressr", "writexl")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)           # resolves file paths relative to the project root (.Rproj)
library(metaselection)  # selection_model(), define_priors() and print() for the fits
library(writexl)        # write_xlsx()
# future and progressr are loaded below, only if refitting

# ---- Data --------------------------------------------------------------------
dat <- readRDS(here("SCData_processed.rds"))
cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- Priors ------------------------------------------------------------------
# Weak default priors from Pustejovsky, Citkowicz, and Joshi (2026), Equations
# (15)-(17). They are nearly flat over the plausible range of each parameter.
# These defaults suit standardized effect sizes such as Fisher's z; for a
# nonstandard effect size (e.g., a raw regression coefficient), the prior on
# beta should be rescaled to the scale of the data.
prior_spec <- define_priors(
  beta_mean = 0, beta_precision = 1 / 16, beta_L = 4,
  tau_mode = 0.2, tau_alpha = 1,
  lambda_mode = 0.5, lambda_precision = 1 / 54, lambda_L = 4
)

# ---- Helper: fit one step-function model ----------------------------------
fit_psm <- function(steps) {
  set.seed(12345)   # same seed as the original script, so results reproduce
  selection_model(
    data = dat, yi = z, sei = sez, cluster = newid,
    selection_type = "step",
    steps          = steps,
    estimator      = "CML",                       # composite marginal likelihood; with the priors above, this is the PML estimator
    priors         = prior_spec,
    bootstrap      = "two-stage",                 # resamples studies, then estimates
    CI_type        = c("large-sample", "percentile"),
    R              = 1999
  )
}

# ---- Fit or load the two models ------------------------------------------------
fit_files <- c("3PSM" = here("sel_3psm_intercept.rds"),
               "4PSM" = here("sel_4psm_intercept.rds"))

if (refit) {
  library(future)     # parallel processing for the bootstrap
  library(progressr)  # progress bar during the bootstrap
  handlers(global = TRUE)
  plan(multisession, workers = 4)   # adjust to the number of cores available

  cat("\nFitting 3PSM (one threshold, p = 0.025)...\n")
  mod_3psm <- fit_psm(0.025)
  cat("\nFitting 4PSM (thresholds p = 0.025 and 0.5)...\n")
  mod_4psm <- fit_psm(c(0.025, 0.5))

  plan(sequential)   # close the parallel workers
  saveRDS(mod_3psm, fit_files["3PSM"])
  saveRDS(mod_4psm, fit_files["4PSM"])
} else {
  missing <- fit_files[!file.exists(fit_files)]
  if (length(missing) > 0)
    stop("Saved fit not found: ", paste(missing, collapse = ", "), ". Set refit <- TRUE.")
  mod_3psm <- readRDS(fit_files["3PSM"])
  mod_4psm <- readRDS(fit_files["4PSM"])
}
print(mod_3psm)
print(mod_4psm)

# ---- Extract estimates ---------------------------------------------------------
# print() with trans_gamma/transf_zeta returns a table with lambda and tau^2
# on their natural scales. Columns: param, Est, SE, CI_lo, CI_hi,
# percentile_lower, percentile_upper.
param_table <- function(mod) {
  res <- NULL
  # capture.output() stops the table being printed again each time it is used
  invisible(capture.output(res <- print(mod, trans_gamma = TRUE, transf_zeta = TRUE)))
  res
}
get_param <- function(res, name) {
  r <- res[res$param == name, ]
  if (nrow(r) == 0) return(NULL)
  r
}

f3 <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "--" else sprintf("%.3f", x)
ci <- function(lo, hi) if (is.null(lo) || length(lo) == 0) "--" else sprintf("[%.3f, %.3f]", lo, hi)

summarise_psm <- function(mod) {
  res <- param_table(mod)
  b  <- get_param(res, "beta")
  t2 <- get_param(res, "tau2")
  l1 <- get_param(res, "lambda1")
  l2 <- get_param(res, "lambda2")   # NULL for 3PSM
  list(
    display = c(
      f3(b$Est),
      ci(b$percentile_lower, b$percentile_upper),
      f3(l1$Est), ci(l1$percentile_lower, l1$percentile_upper),
      if (is.null(l2)) "--" else f3(l2$Est),
      if (is.null(l2)) "--" else ci(l2$percentile_lower, l2$percentile_upper),
      f3(sqrt(t2$Est)),
      nrow(dat), length(unique(dat$newid))
    ),
    values = data.frame(
      Parameter = c("beta", "lambda1", "lambda2"),
      Estimate  = c(b$Est, l1$Est, if (is.null(l2)) NA else l2$Est),
      Lower     = c(b$percentile_lower, l1$percentile_lower,
                    if (is.null(l2)) NA else l2$percentile_lower),
      Upper     = c(b$percentile_upper, l1$percentile_upper,
                    if (is.null(l2)) NA else l2$percentile_upper),
      stringsAsFactors = FALSE
    )
  )
}

s3 <- summarise_psm(mod_3psm)
s4 <- summarise_psm(mod_4psm)

display <- data.frame(
  Parameter = c("Corrected mean", "95% CI (bootstrap)",
                "lambda1", "95% CI (bootstrap)",
                "lambda2", "95% CI (bootstrap)",
                "tau", "Observations", "Studies"),
  `3PSM` = s3$display,
  `4PSM` = s4$display,
  check.names = FALSE, stringsAsFactors = FALSE
)
# lambda1 means different things in 3PSM and 4PSM, so the rows are labelled
# simply lambda1 and lambda2; the note printed below the table in the protocol
# defines each weight for each model.

values <- rbind(cbind(Model = "3PSM", s3$values), cbind(Model = "4PSM", s4$values))

write_xlsx(list(Display = display, Values = values), path = here("Table_SelectionModels.xlsx"))
print(display, row.names = FALSE)
cat("\nSaved Table_SelectionModels.xlsx\n")

elapsed <- proc.time() - start_time
cat(sprintf("Total run time: %.1f minutes.\n", elapsed["elapsed"] / 60))
