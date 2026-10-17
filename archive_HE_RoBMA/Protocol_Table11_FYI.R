library(here)   # resolves file paths relative to the project root (.Rproj)
library(RoBMA)
library(writexl)

start_time <- proc.time()

# ============================================================
# FYI / ROBUSTNESS CHECK -- NOT THE PRIMARY PROTOCOL RESULT
# ============================================================
# Best-practice predictions from the four-moderator FYI model fit in
# Protocol_Table10_FYI.R (Endog_FE, Endog_IV, LaggedDV, Reg_OECDEurope),
# reported alongside Protocol_Table11.R's primary two-moderator predictions
# as a robustness check on how sensitive the results are to moderator
# choice, not as a replacement for them.

# ---- LOAD MODEL ----
# The fitted FYI multilevel RoBMA meta-regression from
# Protocol_Table10_FYI.R. No refitting needed.
cat("Loading Protocol_Table10_FYI.R model...\n")
fit <- readRDS(here("RoBMA_metaregression_FYI.rds"))
cat("Model loaded.\n\n")

# ---- BEST-PRACTICE COVARIATE VALUES ----
# Endog_FE = 1 and Endog_IV = 1 simultaneously: this matches Part 1 Table 7's
# best-practice definition, which follows Xue, Reed, and van Aert (2025) in
# evaluating the idealized combination of both endogeneity-correction
# strategies, even though no single study in the sample uses both at once
# (Endog_IV, Endog_FE, and NoEndogeneity are mutually exclusive in the raw
# data -- see Protocol_DataCleaning.R's Check 4). This is not a data
# inconsistency; it is the same deliberate idealized-profile convention
# already used throughout this protocol's best-practice predictions.
# LaggedDV = 0 matches Part 1 Table 7 as well (no lagged dependent variable).
#
# The two scenarios differ only in region, matching Protocol_Table11.R:
#   BP#1: Non-OECD/Europe (Reg_OECDEurope = 0, implying the combined reference
#         category of Reg_US, Reg_Africa, Reg_Asia, and Reg_Other)
#   BP#2: OECD/Europe     (Reg_OECDEurope = 1)
#
# All binary moderators are encoded as factors with levels c(0, 1) to match
# the factor coding used when the model was fit in Protocol_Table10_FYI.R,
# so that predict() applies the same contrasts.

bp1 <- data.frame(
  Endog_FE       = factor(1, levels = c(0, 1)),
  Endog_IV       = factor(1, levels = c(0, 1)),
  LaggedDV       = factor(0, levels = c(0, 1)),
  Reg_OECDEurope = factor(0, levels = c(0, 1))   # Non-OECD/Europe
)

bp2 <- data.frame(
  Endog_FE       = factor(1, levels = c(0, 1)),
  Endog_IV       = factor(1, levels = c(0, 1)),
  LaggedDV       = factor(0, levels = c(0, 1)),
  Reg_OECDEurope = factor(1, levels = c(0, 1))   # OECD/Europe
)

# ---- PREDICTIONS ----
# Same construction as Protocol_Table11.R: type = "terms" gives the
# credibility interval for the mean (mu, bias-adjusted); type = "estimate"
# adds heterogeneity for the prediction interval.

cat("Computing predictions for BP#1 (Non-OECD/Europe)...\n")
pred_bp1_mean <- predict(fit, newdata = bp1, type = "terms",    bias_adjusted = TRUE)
pred_bp1_pi   <- predict(fit, newdata = bp1, type = "estimate", bias_adjusted = TRUE)

cat("Computing predictions for BP#2 (OECD/Europe)...\n")
pred_bp2_mean <- predict(fit, newdata = bp2, type = "terms",    bias_adjusted = TRUE)
pred_bp2_pi   <- predict(fit, newdata = bp2, type = "estimate", bias_adjusted = TRUE)

# Diagnostic: print summary of predict objects so field names can be verified.
cat("\nSummary of pred_bp1_mean:\n"); print(summary(pred_bp1_mean))
cat("\nSummary of pred_bp1_pi:\n");   print(summary(pred_bp1_pi))

# ---- EXTRACT AND FORMAT RESULTS ----
# Same two-row-per-scenario layout as Protocol_Table11.R.

fmt_est      <- function(x) formatC(x, digits = 3, format = "f")
fmt_interval <- function(lo, hi) paste0("[", fmt_est(lo), ", ", fmt_est(hi), "]")

build_bp_panel <- function(pred_mean, pred_pi) {

  s_mean <- summary(pred_mean)   # type = "terms": posterior for mu
  s_pi   <- summary(pred_pi)     # type = "estimate": posterior for mu + gamma + theta

  z_mean <- s_mean$Mean
  ci_lo  <- s_mean[["0.025"]]
  ci_hi  <- s_mean[["0.975"]]
  pi_lo  <- s_pi[["0.025"]]
  pi_hi  <- s_pi[["0.975"]]

  data.frame(
    check.names = FALSE,
    `Effect Size`      = c("Fisher's z", "PCC"),
    `Mean Prediction`  = c(fmt_est(z_mean), fmt_est(tanh(z_mean))),
    `95% CI`           = c(fmt_interval(ci_lo, ci_hi), fmt_interval(tanh(ci_lo), tanh(ci_hi))),
    `95% PI`           = c(fmt_interval(pi_lo, pi_hi), fmt_interval(tanh(pi_lo), tanh(pi_hi))),
    stringsAsFactors = FALSE
  )
}

cat("\nBuilding BP#1 panel (Non-OECD/Europe)...\n")
panel_bp1 <- build_bp_panel(pred_bp1_mean, pred_bp1_pi)
print(panel_bp1)

cat("\nBuilding BP#2 panel (OECD/Europe)...\n")
panel_bp2 <- build_bp_panel(pred_bp2_mean, pred_bp2_pi)
print(panel_bp2)

# ---- NOTES: "EVALUATED AT" COVARIATE VALUES ----

bp_note <- function(Reg_OECDEurope_bp) {
  paste0(
    "Robustness check, not the primary protocol result. These predictions are derived from ",
    "the four-moderator RoBMA meta-regression in Protocol_Table10_FYI.R (Endog_FE, Endog_IV, LaggedDV, ",
    "Reg_OECDEurope -- the three moderators Table 6 Panel B classifies as \"Evidence of Moderator ",
    "Effect,\" plus Reg_OECDEurope for this regional comparison), rather than the protocol's primary ",
    "BIC-selected two-moderator model (Table 10 / Table 11). Predictions are evaluated at Endog_FE = 1, ",
    "Endog_IV = 1, LaggedDV = 0, Reg_OECDEurope = ", Reg_OECDEurope_bp, " ",
    "(the reference category -- combining Reg_US, Reg_Africa, Reg_Asia, and Reg_Other -- is implied ",
    "when Reg_OECDEurope = 0). se(z) and publication year are not covariates in this model: RoBMA ",
    "corrects for publication bias internally rather than through a linear se(z) term, and publication ",
    "year was not retained by the BIC pre-selection in Table 9."
  )
}

note_bp1 <- bp_note(Reg_OECDEurope_bp = 0)
note_bp2 <- bp_note(Reg_OECDEurope_bp = 1)

cat("\nBP#1 note:\n", note_bp1, "\n")
cat("\nBP#2 note:\n", note_bp2, "\n")

notes_sheet <- data.frame(
  Scenario = c("BP1 - Non-OECD/Europe", "BP2 - OECD/Europe"),
  Note     = c(note_bp1, note_bp2),
  stringsAsFactors = FALSE
)

# ---- EXPORT TO EXCEL ----
# Same three-sheet structure as Table11_RoBMA_BestPractice.xlsx, saved under
# a distinct filename so it does not overwrite the primary protocol result.

write_xlsx(
  list(
    "BP1 - Non-OECD Europe" = panel_bp1,
    "BP2 - OECD Europe"     = panel_bp2,
    "Notes"                 = notes_sheet
  ),
  here("Table11_FYI_BestPractice.xlsx")
)
cat("\nResults saved to Table11_FYI_BestPractice.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds\n", total_time["elapsed"]))
cat("\nProtocol_Table11_FYI.R complete.\n")
