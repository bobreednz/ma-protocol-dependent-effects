# =============================================================================
# Protocol_AppendixPubBias.R
# =============================================================================
# Purpose: Robustness checks (Section 6) -- alternative publication bias
# estimators.
#
# Produces a single table with four columns, all estimating the average effect
# size adjusted for publication bias. This is the robustness-check complement
# to Table 5 (FAT-PET under standard CHE weighting).
#
#   Column 1: CHE-ISCW-PET-PEESE, no controls
#             PET-PEESE under inverse sampling covariance weighting, with
#             only the bias proxy as a regressor (analogous to Table 5
#             Panel A).
#
#   Column 2: CHE-ISCW-PET-PEESE, with controls
#             Same as column 1 but with all 18 study characteristics added
#             as regressors, centered at their sample means (analogous to
#             Table 5 Panel B). The intercept is the "Effect beyond bias"
#             at mean covariate values.
#
#   Column 3: PML one-step selection model (3PSM)
#             Single threshold at one-sided p = 0.025. Estimates the
#             relative selection probability lambda1 for non-significant
#             results versus significant ones.
#
#   Column 4: PML two-step selection model (4PSM)
#             Two thresholds at one-sided p = 0.025 and 0.5. Adds lambda2
#             for results in the opposite direction.
#
# PET-PEESE decision rule (applied separately for each CHE-ISCW column):
#   Use PEESE (sez^2) if the PET intercept is significant at p < 0.10 AND
#   the PET intercept is positive. Otherwise use PET (sez). A row in the
#   table records which was selected.
#
# CHE-ISCW weights by the inverse of the sampling covariance matrix rather
# than the inverse of the total (sampling + random effects) covariance matrix.
# This makes the estimator less susceptible to publication bias because the
# weights do not depend on variance components that are inflated by selective
# reporting (Chen & Pustejovsky, 2026).
#
# PML selection models use two-stage cluster bootstrap (R = 1999) for CIs,
# which respects the dependent data structure. They are estimated without
# controls because adding many moderators causes convergence instability in
# bootstrap iterations and because the intercept-only selection model is
# standard in the literature (Pustejovsky, Citkowicz, & Joshi, 2026).
#
# Outputs:
#   TableA_PubBias.xlsx
#   FigureA_SelectionPlot_3PSM.png
#   FigureA_SelectionPlot_4PSM.png
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(here)            # resolves file paths relative to the project root (.Rproj)
library(metafor)        # rma.mv() for CHE-ISCW models
library(clubSandwich)   # robust() for CR2 standard errors
library(metaselection)  # selection_model() and define_priors()
library(future)         # plan(multisession) for parallel bootstrapping
library(progressr)      # progress bar during bootstrap
library(tidyverse)      # data manipulation
library(writexl)        # write_xlsx()

# If metaselection is not yet installed (GitHub only, not on CRAN):
# remotes::install_github("jepusto/metaselection")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ── Data ──────────────────────────────────────────────────────────────────────

dat <- readRDS(here("SCData_processed.rds"))

cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ── Sampling covariance matrix and ISCW weight matrix ─────────────────────────

# Imputed block-diagonal covariance matrix with rho = 0.5 (as used throughout
# the protocol). metafor::vcalc() replaced the deprecated
# clubSandwich::impute_covariance_matrix().
V_mat <- vcalc(
  vi      = sez^2,
  cluster = newid,
  obs     = obs,
  rho     = 0.5,
  data    = dat
)

# The ISCW weight matrix is the inverse of V_mat. Passing this as W to
# rma.mv() makes the estimator weight by sampling precision alone rather than
# total precision (sampling + random effects), reducing susceptibility to
# publication bias. vcalc() returns a single block-diagonal matrix, so it can
# be inverted directly.
W_iscw <- solve(V_mat)

# ── Center the 18 control variables ───────────────────────────────────────────

# Centering at sample means means the intercept in every with-controls model
# equals the predicted effect when sez = 0 (no bias) and all controls are at
# their sample means. This is the same centering strategy used in Table 5
# Panel B. Continuous variables (PubYear, NumberSCVars) are centered at their
# means; binary variables are also centered, which is equivalent to computing
# the average marginal prediction.

controls <- c("SC1_Cognitive", "SC1_Structural", "DV_GrowthRate",
              "PubYear", "Published", "LaggedDV", "LaggedSC",
              "NumberSCVars", "Endog_IV", "Endog_FE",
              "CityLevel", "RegionLevel", "CountryLevel",
              "PanelData", "Reg_OECDEurope", "Reg_US",
              "Reg_Africa", "Reg_Asia")

for (v in controls) {
  dat[[paste0(v, "_c")]] <- dat[[v]] - mean(dat[[v]], na.rm = TRUE)
}

# Build a formula string for the 18 centered controls (used in with-controls
# models). We paste them together so we can drop them into a formula easily.
controls_rhs <- paste(paste0(controls, "_c"), collapse = " + ")

# =============================================================================
# Part 1: CHE-ISCW-PET-PEESE (no controls)
# =============================================================================

cat("\nFitting CHE-ISCW-PET (no controls)...\n")

# PET: bias proxy is the standard error (sez). The intercept gives the
# estimated mean effect when sez = 0, i.e., in a hypothetical study with no
# sampling error -- the "Effect beyond bias."
pet_nc <- rma.mv(
  yi     = z,
  V      = V_mat,
  mods   = ~ sez,            # sez is the bias proxy; not centered
  W      = W_iscw,           # ISCW weight matrix
  random = ~ 1 | newid/obs,
  data   = dat,
  method = "REML"
)
pet_nc_rob <- robust(pet_nc, cluster = dat$newid, clubSandwich = TRUE)

# PEESE: bias proxy is sez^2. Use I() to square within the formula.
cat("Fitting CHE-ISCW-PEESE (no controls)...\n")
peese_nc <- rma.mv(
  yi     = z,
  V      = V_mat,
  mods   = ~ I(sez^2),
  W      = W_iscw,
  random = ~ 1 | newid/obs,
  data   = dat,
  method = "REML"
)
peese_nc_rob <- robust(peese_nc, cluster = dat$newid, clubSandwich = TRUE)

# PET-PEESE decision rule: use PEESE if PET intercept is significant at
# p < 0.10 AND positive; otherwise report PET.
# pval[1] is the two-sided p-value for the intercept from robust().
use_peese_nc <- (pet_nc_rob$pval[1] < 0.10) & (pet_nc_rob$beta[1] > 0)
selected_nc     <- if (use_peese_nc) peese_nc_rob else pet_nc_rob
proxy_label_nc  <- if (use_peese_nc) "PEESE (sez^2)" else "PET (sez)"

cat(sprintf("No-controls: %s selected (PET intercept p = %.3f, beta = %.4f).\n",
            proxy_label_nc, pet_nc_rob$pval[1], pet_nc_rob$beta[1]))

# =============================================================================
# Part 2: CHE-ISCW-PET-PEESE (with controls)
# =============================================================================

cat("\nFitting CHE-ISCW-PET (with controls)...\n")

pet_wc_formula   <- as.formula(paste("~ sez +", controls_rhs))
peese_wc_formula <- as.formula(paste("~ I(sez^2) +", controls_rhs))

pet_wc <- rma.mv(
  yi     = z,
  V      = V_mat,
  mods   = pet_wc_formula,
  W      = W_iscw,
  random = ~ 1 | newid/obs,
  data   = dat,
  method = "REML"
)
pet_wc_rob <- robust(pet_wc, cluster = dat$newid, clubSandwich = TRUE)

cat("Fitting CHE-ISCW-PEESE (with controls)...\n")
peese_wc <- rma.mv(
  yi     = z,
  V      = V_mat,
  mods   = peese_wc_formula,
  W      = W_iscw,
  random = ~ 1 | newid/obs,
  data   = dat,
  method = "REML"
)
peese_wc_rob <- robust(peese_wc, cluster = dat$newid, clubSandwich = TRUE)

# Decision rule: applied separately for the with-controls version.
use_peese_wc   <- (pet_wc_rob$pval[1] < 0.10) & (pet_wc_rob$beta[1] > 0)
selected_wc    <- if (use_peese_wc) peese_wc_rob else pet_wc_rob
proxy_label_wc <- if (use_peese_wc) "PEESE (sez^2)" else "PET (sez)"

cat(sprintf("With-controls: %s selected (PET intercept p = %.3f, beta = %.4f).\n",
            proxy_label_wc, pet_wc_rob$pval[1], pet_wc_rob$beta[1]))

# =============================================================================
# Part 3: Prior specification for PML selection models
# =============================================================================

# Weak default priors from Pustejovsky et al. (2026), Equations (15)-(17).
# Nearly flat over the plausible parameter range; only informative at extremes.
prior_spec <- define_priors(
  beta_mean        = 0,
  beta_precision   = 1 / 16,
  beta_L           = 4,
  tau_mode         = 0.2,
  tau_alpha        = 1,
  lambda_mode      = 0.5,
  lambda_precision = 1 / 54,
  lambda_L         = 4
)

# Enable parallel processing and a progress bar for the bootstrap loop.
# Adjust workers to match the number of cores available on your machine.
handlers(global = TRUE)
plan(multisession, workers = 4)

# =============================================================================
# Part 4: One-step PML selection model (3PSM)
# =============================================================================
#
# Single threshold at one-sided p = 0.025 (two-sided p = 0.05 for positive
# effects). lambda1 is the relative selection probability for non-significant
# or negative results versus statistically significant positive results.
# lambda1 = 1 means no selection; lambda1 < 1 means non-significant results
# are under-represented.

cat("\nFitting one-step PML selection model (3PSM)...\n")
cat("(Bootstrap R = 1999; allow 15-30 minutes.)\n")

set.seed(12345)
mod_3psm <- selection_model(
  data           = dat,
  yi             = z,
  sei            = sez,
  cluster        = newid,
  selection_type = "step",
  steps          = 0.025,
  estimator      = "CML",       # CML = PML (penalized marginal likelihood)
  priors         = prior_spec,
  bootstrap      = "two-stage", # respects study-level clustering
  CI_type        = c("large-sample", "percentile"),
  R              = 1999
)

print(mod_3psm)

# Save the selection function plot. transform = "sqrt" spreads the lower
# end of the lambda range for readability.
png(here("FigureA_SelectionPlot_3PSM.png"), width = 6, height = 5, units = "in", res = 300)
print(selection_plot(mod_3psm, transform = "sqrt"))
dev.off()
cat("FigureA_SelectionPlot_3PSM.png saved.\n")

# =============================================================================
# Part 5: Two-step PML selection model (4PSM)
# =============================================================================
#
# Two thresholds distinguish three regions of the one-sided p-value:
#   p <= 0.025       : significant positive results (reference, lambda = 1)
#   0.025 < p <= 0.5 : positive but not significant (lambda1)
#   p > 0.5          : negative results (lambda2)
# lambda2 is often estimated imprecisely when negative results are rare.

cat("\nFitting two-step PML selection model (4PSM)...\n")
cat("(Bootstrap R = 1999; allow 15-30 minutes.)\n")

set.seed(12345)
mod_4psm <- selection_model(
  data           = dat,
  yi             = z,
  sei            = sez,
  cluster        = newid,
  selection_type = "step",
  steps          = c(0.025, 0.5),
  estimator      = "CML",
  priors         = prior_spec,
  bootstrap      = "two-stage",
  CI_type        = c("large-sample", "percentile"),
  R              = 1999
)

print(mod_4psm)

png(here("FigureA_SelectionPlot_4PSM.png"), width = 6, height = 5, units = "in", res = 300)
print(selection_plot(mod_4psm, transform = "sqrt"))
dev.off()
cat("FigureA_SelectionPlot_4PSM.png saved.\n")

# Close parallel workers now that bootstrapping is done.
plan(sequential)

# =============================================================================
# Part 6: Assemble the results table
# =============================================================================

# ── 6a. Helper: extract intercept and bias coefficient from robust() output ───

# robust() objects from metafor have: beta (matrix), se (vector), ci.lb,
# ci.ub (Satterthwaite-based CIs), pval, dfs. Row 1 is the intercept ("intrcpt");
# row 2 is the bias proxy coefficient (sez or sez^2).

extract_iscw <- function(rob_mod, bias_label, tau, omega) {
  list(
    effect     = round(as.numeric(rob_mod$beta[1]), 4),
    se         = round(rob_mod$se[1], 4),
    ci         = paste0("[", round(rob_mod$ci.lb[1], 4), ", ",
                             round(rob_mod$ci.ub[1], 4), "]"),
    proxy      = bias_label,
    bias_coef  = round(as.numeric(rob_mod$beta[2]), 4),
    bias_se    = round(rob_mod$se[2], 4),
    tau        = round(sqrt(tau), 4),
    omega      = round(sqrt(omega), 4)
  )
}

# Retrieve sigma2 from the underlying rma.mv object (not the robust wrapper).
# For the no-controls models:
sigma2_nc <- if (use_peese_nc) peese_nc$sigma2 else pet_nc$sigma2
sigma2_wc <- if (use_peese_wc) peese_wc$sigma2 else pet_wc$sigma2

res_nc <- extract_iscw(selected_nc, proxy_label_nc,
                       tau   = sigma2_nc[1],
                       omega = sigma2_nc[2])
res_wc <- extract_iscw(selected_wc, proxy_label_wc,
                       tau   = sigma2_wc[1],
                       omega = sigma2_wc[2])

# ── 6b. Helper: extract results from PML selection model print output ─────────

# print() with trans_gamma = TRUE converts log(tau^2) -> tau^2.
# transf_zeta = TRUE converts log(lambda) -> lambda.
# Columns: param, Est, SE, CI_lo, CI_hi, percentile_lower, percentile_upper.

extract_pml <- function(mod) {
  res <- print(mod, trans_gamma = TRUE, transf_zeta = TRUE)

  get_row <- function(param_name) {
    r <- res[res$param == param_name, ]
    if (nrow(r) == 0) return(NULL)
    r
  }

  beta_row    <- get_row("beta")
  tau2_row    <- get_row("tau2")
  lambda1_row <- get_row("lambda1")
  lambda2_row <- get_row("lambda2")

  fmt_ci <- function(lo, hi) paste0("[", round(lo, 4), ", ", round(hi, 4), "]")

  list(
    effect          = round(beta_row$Est, 4),
    se              = round(beta_row$SE, 4),
    ci_ls           = fmt_ci(beta_row$CI_lo, beta_row$CI_hi),
    ci_boot         = fmt_ci(beta_row$percentile_lower, beta_row$percentile_upper),
    lambda1         = round(lambda1_row$Est, 4),
    lambda1_se      = round(lambda1_row$SE, 4),
    lambda1_ci_ls   = fmt_ci(lambda1_row$CI_lo, lambda1_row$CI_hi),
    lambda1_ci_boot = fmt_ci(lambda1_row$percentile_lower, lambda1_row$percentile_upper),
    lambda2         = if (!is.null(lambda2_row)) round(lambda2_row$Est, 4)          else NA,
    lambda2_se      = if (!is.null(lambda2_row)) round(lambda2_row$SE, 4)           else NA,
    lambda2_ci_ls   = if (!is.null(lambda2_row)) fmt_ci(lambda2_row$CI_lo,
                                                         lambda2_row$CI_hi)         else NA,
    lambda2_ci_boot = if (!is.null(lambda2_row)) fmt_ci(lambda2_row$percentile_lower,
                                                         lambda2_row$percentile_upper) else NA,
    tau             = round(sqrt(tau2_row$Est), 4)
  )
}

res_3psm <- extract_pml(mod_3psm)
res_4psm <- extract_pml(mod_4psm)

# ── 6c. Build the table ───────────────────────────────────────────────────────

# Rows common to all four columns are stacked. Rows specific to CHE-ISCW
# (bias coefficient, omega) and to PML (lambda parameters, bootstrap CIs)
# are filled with NA where not applicable. A "PET or PEESE" row records
# which bias proxy was selected for each CHE-ISCW column.

tbl <- data.frame(
  Parameter = c(
    # Effect estimate block
    "Effect beyond bias (z)",
    "SE",
    "95% CI (Satterthwaite)",        # for CHE-ISCW columns
    "95% CI (large-sample)",         # for PML columns
    "95% CI (percentile bootstrap)", # for PML columns
    # Bias proxy block
    "Bias proxy selected",
    "Bias coefficient",
    "SE (bias coefficient)",
    # Lambda block (PML only)
    "lambda1",
    "SE (lambda1)",
    "95% CI lambda1 (large-sample)",
    "95% CI lambda1 (bootstrap)",
    "lambda2",
    "SE (lambda2)",
    "95% CI lambda2 (large-sample)",
    "95% CI lambda2 (bootstrap)",
    # Variance components
    "tau",
    "omega",
    # Sample info
    "Observations",
    "Studies"
  ),

  # Column 1: CHE-ISCW-PET-PEESE, no controls
  CHE_ISCW_no_controls = c(
    res_nc$effect, res_nc$se, res_nc$ci,
    NA, NA,                          # no large-sample / bootstrap CI for CHE-ISCW
    res_nc$proxy,
    res_nc$bias_coef, res_nc$bias_se,
    NA, NA, NA, NA,                  # no lambda for CHE-ISCW
    NA, NA, NA, NA,
    res_nc$tau, res_nc$omega,
    nrow(dat), length(unique(dat$newid))
  ),

  # Column 2: CHE-ISCW-PET-PEESE, with controls
  CHE_ISCW_with_controls = c(
    res_wc$effect, res_wc$se, res_wc$ci,
    NA, NA,
    res_wc$proxy,
    res_wc$bias_coef, res_wc$bias_se,
    NA, NA, NA, NA,
    NA, NA, NA, NA,
    res_wc$tau, res_wc$omega,
    nrow(dat), length(unique(dat$newid))
  ),

  # Column 3: PML one-step (3PSM)
  PML_1step = c(
    res_3psm$effect, res_3psm$se,
    NA,                              # no Satterthwaite CI for PML
    res_3psm$ci_ls, res_3psm$ci_boot,
    NA, NA, NA,                      # no bias proxy row for PML
    res_3psm$lambda1, res_3psm$lambda1_se,
    res_3psm$lambda1_ci_ls, res_3psm$lambda1_ci_boot,
    NA, NA, NA, NA,                  # no lambda2 for 3PSM
    res_3psm$tau,
    NA,                              # no omega in marginal selection model
    nrow(dat), length(unique(dat$newid))
  ),

  # Column 4: PML two-step (4PSM)
  PML_2step = c(
    res_4psm$effect, res_4psm$se,
    NA,
    res_4psm$ci_ls, res_4psm$ci_boot,
    NA, NA, NA,
    res_4psm$lambda1, res_4psm$lambda1_se,
    res_4psm$lambda1_ci_ls, res_4psm$lambda1_ci_boot,
    res_4psm$lambda2, res_4psm$lambda2_se,
    res_4psm$lambda2_ci_ls, res_4psm$lambda2_ci_boot,
    res_4psm$tau,
    NA,
    nrow(dat), length(unique(dat$newid))
  ),

  stringsAsFactors = FALSE
)

# =============================================================================
# Part 7: Save outputs
# =============================================================================

write_xlsx(
  list("Alt Pub Bias" = tbl),
  path = here("TableA_PubBias.xlsx")
)

cat("\nAll outputs saved:\n")
cat("  TableA_PubBias.xlsx\n")
cat("  FigureA_SelectionPlot_3PSM.png\n")
cat("  FigureA_SelectionPlot_4PSM.png\n")
cat(sprintf("\nPET-PEESE selection summary:\n"))
cat(sprintf("  No-controls:   %s\n", proxy_label_nc))
cat(sprintf("  With-controls: %s\n", proxy_label_wc))

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
