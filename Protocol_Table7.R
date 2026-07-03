setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("metafor", "clubSandwich", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(metafor)       # meta-analytic models via rma.mv()
library(clubSandwich)  # CR2 clustered standard errors
library(tidyverse)     # data manipulation
library(writexl)       # export to Excel

# ============================================================
# Load data
# ============================================================

DT <- readRDS("SCData_processed.rds")

# ============================================================
# Impute block-diagonal covariance matrix (for CHE)
# ============================================================

rho   <- 0.5

V_mat <- vcalc(
  vi      = sez^2,      # sampling variances
  cluster = newid,      # study identifier
  obs     = obs,        # observation identifier within study
  rho     = rho,        # assumed within-study correlation
  data    = DT
)

# ============================================================
# Define best practice values
#
# "Best practice" means methodologically ideal study
# characteristics. Following Table 11 of Xue, Reed, and van
# Aert (2024), the best practice study:
#
#   - Reports effects for economic growth rates (not GDP levels)
#   - Uses cognitive social capital
#   - Has no publication bias (se(z) = 0)
#   - Is a published article
#   - Does not use lagged DV or lagged SC
#   - Includes only one social capital variable
#   - Uses both instrumental variables and fixed effects
#     to address endogeneity
#   - Uses country-level panel data
#   - Is published in the sample's mean publication year
#
# The two scenarios differ only in region:
#   BP#1: Non-OECD/Europe (all region dummies = 0)
#   BP#2: OECD/Europe    (Reg_OECDEurope = 1)
#
# PubYear is set to the sample mean (continuous variable).
# NumberSCVars is set to 1 (one social capital variable).
# ============================================================

bp_PubYear      <- mean(DT$PubYear)   # continuous: use sample mean
bp_NumberSCVars <- 1                  # single social capital variable

# ============================================================
# Center variables at best practice values
#
# This is the key step. Subtracting a variable's best practice
# value before estimation means the intercept in any model
# equals the predicted effect when all centered variables are
# zero -- i.e., exactly at best practice. The CR2 standard
# error and Satterthwaite df on the intercept then directly
# characterize uncertainty about the best practice estimate,
# with no post-estimation prediction or delta-method required.
#
# We define a small helper function that does the centering
# for a given value of Reg_OECDEurope (the only variable that
# differs between the two scenarios).
# ============================================================

center_at_bp <- function(data, Reg_OECDEurope_bp) {
  data %>%
    mutate(
      # se(z) = 0 at best practice, so sez_c = sez - 0 = sez.
      # The centered name is kept for clarity in the formula.
      sez_c            = sez            - 0,

      # Social capital type: cognitive = 1, structural = 0
      SC1_Cognitive_c  = SC1_Cognitive  - 1,
      SC1_Structural_c = SC1_Structural - 0,

      # Outcome: growth rate = 1
      DV_GrowthRate_c  = DV_GrowthRate  - 1,

      # Publication year: sample mean (continuous variable)
      PubYear_c        = PubYear        - bp_PubYear,

      # Published study
      Published_c      = Published      - 1,

      # No lagged DV or SC
      LaggedDV_c       = LaggedDV       - 0,
      LaggedSC_c       = LaggedSC       - 0,

      # Single social capital variable (continuous)
      NumberSCVars_c   = NumberSCVars   - bp_NumberSCVars,

      # IV and FE used for endogeneity
      Endog_IV_c       = Endog_IV       - 1,
      Endog_FE_c       = Endog_FE       - 1,

      # Country-level data (not city or region level)
      CityLevel_c      = CityLevel      - 0,
      RegionLevel_c    = RegionLevel    - 0,
      CountryLevel_c   = CountryLevel   - 1,

      # Panel data
      PanelData_c      = PanelData      - 1,

      # Region: the one variable that differs between scenarios
      Reg_OECDEurope_c = Reg_OECDEurope - Reg_OECDEurope_bp,

      # Other regions: 0 at best practice
      Reg_US_c         = Reg_US         - 0,
      Reg_Africa_c     = Reg_Africa     - 0,
      Reg_Asia_c       = Reg_Asia       - 0
    )
}

# Create two centered datasets
DT_bp1 <- center_at_bp(DT, Reg_OECDEurope_bp = 0)   # BP#1: non-OECD/Europe
DT_bp2 <- center_at_bp(DT, Reg_OECDEurope_bp = 1)   # BP#2: OECD/Europe

# ============================================================
# Model formula
#
# Same 19 regressors as the CHE model in Table 6 (Panel A) and Table 10, but now
# using the centered variables. Because sez_c = sez (best
# practice sez = 0), and centering does not change sampling
# variances, V = sez^2 still holds.
# ============================================================

formula_bp <- z ~ 1 + sez_c +
  SC1_Cognitive_c + SC1_Structural_c + DV_GrowthRate_c +
  PubYear_c + Published_c + LaggedDV_c + LaggedSC_c + NumberSCVars_c +
  Endog_IV_c + Endog_FE_c + CityLevel_c + RegionLevel_c + CountryLevel_c +
  PanelData_c + Reg_OECDEurope_c + Reg_US_c + Reg_Africa_c + Reg_Asia_c

# ============================================================
# Helper functions
# ============================================================

# Significance stars
stars <- function(p) {
  case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ "")
}

# Format an estimate with significance stars
fmt_est <- function(est, p) paste0(formatC(est, digits = 3, format = "f"), stars(p))

# Format a confidence or prediction interval as [lo, hi]
fmt_interval <- function(lo, hi) {
  paste0("[",
         formatC(lo, digits = 3, format = "f"),
         ", ",
         formatC(hi, digits = 3, format = "f"),
         "]")
}

# ============================================================
# Core function: fit the CHE model for one BP scenario and
# return a formatted table
#
# Only CHE is fit here. CHE was established as the preferred
# estimator in the Section 5 Decision Point and as the sole
# inferential model for the meta-regression in Table 6 ("BMA
# informs interpretation; CHE determines estimation and
# inference"). Table 7 applies that same chosen model to the
# best-practice covariate values, so there is no remaining
# role for FE/RE/CE/HE here -- unlike Table 4, this is not an
# estimator-comparison step.
#
# The returned table has rows for scale (Fisher's z, PCC) and
# columns for Mean Prediction, 95% CI, and 95% PI -- a two-row
# panel rather than the old five-model-wide layout, since a
# single-model table has nothing left to compare across columns.
#
# Arguments:
#   DT_bp  -- dataset with all variables centered at BP values
#
# Steps:
#   1. Fit CHE via REML
#   2. Compute CR2 standard errors via coef_test()
#   3. Extract the intercept row (= BP estimate and its SE)
#   4. Compute 95% confidence interval (CI) using the
#      Satterthwaite df from coef_test()
#   5. Compute 95% prediction interval (PI) by adding the
#      model's variance components to the CI standard error.
#      The PI answers: what range of effects would we expect
#      from a new best-practice study? It is wider than the
#      CI because it includes residual between-study
#      heterogeneity (tau^2 and omega^2).
#   6. Convert the Fisher's z estimate and both intervals to PCC
#      via tanh(); tanh() is monotonic, so applying it directly
#      to the z-scale interval endpoints gives the correct PCC
#      interval.
# ============================================================

build_bp_panel <- function(DT_bp) {

  # --- Fit CHE via REML ---
  # V_mat is the global imputed block-diagonal matrix.
  # Centering does not change sez, so V = sez^2 is unchanged.

  m_CHE <- rma.mv(formula_bp,
                   V      = V_mat,
                   random = ~ 1 | newid/obs,
                   data   = DT_bp,
                   method = "REML")

  # --- CR2 standard errors, clustered at the study level ---

  ct_CHE <- coef_test(m_CHE, cluster = DT_bp$newid, vcov = "CR2")

  # --- Extract the intercept row ---
  # The intercept is the best practice estimate because all
  # centered variables equal zero at best practice values.

  int_CHE <- ct_CHE[ct_CHE$Coef == "intrcpt", ]

  # --- Helper: compute CI bounds ---
  # Uses the Satterthwaite df from coef_test() (column df_Satt).
  # This df is specific to the intercept coefficient and accounts
  # for the cluster-robust variance estimation.

  ci_bounds <- function(int_row) {
    t_crit <- qt(0.975, df = int_row$df_Satt)
    c(lo = int_row$beta - t_crit * int_row$SE,
      hi = int_row$beta + t_crit * int_row$SE)
  }

  # --- Helper: compute PI bounds ---
  # The prediction interval adds residual heterogeneity to the
  # SE before computing the interval. sigma2_total is the sum
  # of all variance components from the fitted model.

  pi_bounds <- function(int_row, sigma2_total) {
    t_crit  <- qt(0.975, df = int_row$df_Satt)
    pred_se <- sqrt(int_row$SE^2 + sigma2_total)
    c(lo = int_row$beta - t_crit * pred_se,
      hi = int_row$beta + t_crit * pred_se)
  }

  # Compute the CI for CHE

  ci_CHE <- ci_bounds(int_CHE)

  # Compute the PI for CHE
  # CHE has two variance components (tau^2 for study, omega^2 for
  # observation within study), so sigma2_total is their sum.

  pi_CHE <- pi_bounds(int_CHE, sigma2_total = sum(m_CHE$sigma2))

  # --- Assemble and return the formatted table ---
  # Row 1 is the Fisher's z scale (the scale the model is fit on);
  # row 2 converts the same estimate and intervals to PCC via
  # tanh(). CI reflects estimation uncertainty about the mean;
  # PI additionally incorporates residual heterogeneity.

  tibble(
    `Effect Size` = c("Fisher's z", "PCC"),

    `Mean Prediction` = c(
      fmt_est(int_CHE$beta, int_CHE$p_Satt),
      formatC(tanh(int_CHE$beta), digits = 3, format = "f")
    ),

    `95% CI` = c(
      fmt_interval(ci_CHE["lo"], ci_CHE["hi"]),
      fmt_interval(tanh(ci_CHE["lo"]), tanh(ci_CHE["hi"]))
    ),

    `95% PI` = c(
      fmt_interval(pi_CHE["lo"], pi_CHE["hi"]),
      fmt_interval(tanh(pi_CHE["lo"]), tanh(pi_CHE["hi"]))
    )
  )
}

# ============================================================
# Build tables for both best practice scenarios
#
# BP#1 (non-OECD/Europe) is listed first throughout the protocol
# (Table 7 and its RoBMA counterpart, Table 11), so that order is
# preserved here.
# ============================================================

cat("Fitting CHE model for Best Practice #1 (non-OECD/Europe)...\n")
panel_bp1 <- build_bp_panel(DT_bp1)

cat("Fitting CHE model for Best Practice #2 (OECD/Europe)...\n")
panel_bp2 <- build_bp_panel(DT_bp2)

# ============================================================
# Build a complete "evaluated at" note for each scenario
#
# Earlier versions of this table appended all 19 best-practice
# covariate values as extra rows below the estimates. With the
# table reduced to CHE only, that full list is now written as a
# single note per scenario instead of extra rows, but every value
# is still listed explicitly (no truncation) so the table remains
# fully self-documenting. The notes are printed to console and
# saved to a "Notes" sheet, following the same console-print
# convention used for Table 3's footnotes.
# ============================================================

bp_note <- function(Reg_OECDEurope_bp) {
  paste0(
    "\"Best Practice\" predictions are derived from the CHE estimates in Table 6, Panel A. ",
    "Predictions are evaluated at se(z) = 0, SC1_Cognitive = 1, SC1_Structural = 0, ",
    "SC1_Other = 0 (reference), DV_GrowthRate = 1, DV_GDPLevel = 0 (reference), PubYear = ",
    round(bp_PubYear, 1), " (sample mean), Published = 1, LaggedDV = 0, LaggedSC = 0, NumberSCVars = ",
    bp_NumberSCVars, ", Endog_IV = 1, Endog_FE = 1, CityLevel = 0, RegionLevel = 0, CountryLevel = 1, ",
    "PanelData = 1, Reg_OECDEurope = ", Reg_OECDEurope_bp, ", Reg_US = 0, Reg_Africa = 0, Reg_Asia = 0."
  )
}

note_bp1 <- bp_note(Reg_OECDEurope_bp = 0)
note_bp2 <- bp_note(Reg_OECDEurope_bp = 1)

cat("\nBP#1 note:\n", note_bp1, "\n")
cat("\nBP#2 note:\n", note_bp2, "\n")

notes_sheet <- tibble(
  Scenario = c("BP1 - Non-OECD/Europe", "BP2 - OECD/Europe"),
  Note     = c(note_bp1, note_bp2)
)

# ============================================================
# Export both panels and the notes to Excel
# ============================================================

write_xlsx(
  list(
    "BP1 - Non-OECD Europe" = panel_bp1,
    "BP2 - OECD Europe"     = panel_bp2,
    "Notes"                 = notes_sheet
  ),
  "Table7_Protocol.xlsx"
)

cat("Done. Output saved to Table7_Protocol.xlsx\n")

# -- Run time ------------------------------------------------------------------
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
