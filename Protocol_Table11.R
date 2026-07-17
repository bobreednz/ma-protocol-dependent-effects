library(RoBMA)
library(writexl)

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

start_time <- proc.time()

# ---- LOAD MODEL ----
# The fitted multilevel RoBMA meta-regression from Protocol_Table10.R. No refitting needed.
cat("Loading Protocol_Table10.R model...\n")
fit <- readRDS("RoBMA_metaregression.rds")
cat("Model loaded.\n\n")

# ---- BEST-PRACTICE COVARIATE VALUES ----
# Following Part 1 (Table 7) and Xue, Reed, and van Aert (2025), the
# best-practice study uses cognitive social capital, no lagged DV, one social
# capital variable, and fixed effects for endogeneity. sez and PubYear are not
# in the model (RoBMA corrects for publication bias internally; PubYear was not
# retained by Table 9's BIC pre-selection).
#
# The two scenarios differ only in region:
#   BP#1: Non-OECD/Europe (all four region dummies = 0, implying Reg_Other)
#   BP#2: OECD/Europe     (Reg_OECDEurope = 1, all other region dummies = 0)
#
# NumberSCVars is passed as raw value (1); RoBMA standardizes internally.

bp1 <- data.frame(
  SC1_Cognitive  = 1,
  LaggedDV       = 0,
  NumberSCVars   = 1,
  Endog_FE       = 1,
  Reg_OECDEurope = 0,   # Non-OECD/Europe: all four region dummies = 0
  Reg_US         = 0,
  Reg_Africa     = 0,
  Reg_Asia       = 0
)

bp2 <- data.frame(
  SC1_Cognitive  = 1,
  LaggedDV       = 0,
  NumberSCVars   = 1,
  Endog_FE       = 1,
  Reg_OECDEurope = 1,   # OECD/Europe
  Reg_US         = 0,
  Reg_Africa     = 0,
  Reg_Asia       = 0
)

# ---- PREDICTIONS ----
# predict() with type = "terms" returns the model-averaged posterior for the
# mean effect (mu) at the specified covariate values, corrected for publication
# bias (bias_adjusted = TRUE removes PET/PEESE terms from mu). This is
# analogous to the confidence interval row in Part 1 Table 7: it captures
# uncertainty about the mean effect, not variability across future studies.
#
# type = "estimate" adds random-effect heterogeneity (mu + gamma + theta),
# giving a prediction interval: the range of true effects we would expect
# from a new best-practice study. This is analogous to the PI row in Part 1.
#
# conditional = FALSE (default) gives fully model-averaged predictions.

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
# Following the same two-row-per-scenario layout as Part 1 Table 7 (see
# build_bp_panel() in Protocol_Table7.R): rows are Fisher's z and PCC;
# columns are Mean Prediction, 95% CI, and 95% PI. The mean prediction is
# taken from the CrI-for-mu prediction (type = "terms") and shared across
# both interval columns, since mu's posterior mean is essentially
# unchanged by adding heterogeneity; only the interval width differs
# between the CrI (mu only) and the PI (mu + gamma + theta).

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
# Following the same per-scenario note convention as Part 1 Table 7
# (bp_note() in Protocol_Table7.R): every covariate value is listed
# explicitly rather than truncated, and the note is printed to console
# and saved to its own "Notes" sheet rather than appended as table rows.
#
# RoBMA's meta-regression (Table 10) uses only the 8 moderators retained
# after BIC pre-selection (Table 9) plus the four region dummies. Unlike
# Part 1 Table 7, se(z) and publication year are not covariates here:
# RoBMA corrects for publication bias internally rather than through a
# linear se(z) term, and publication year was not among the moderators
# retained by Table 9's BIC pre-selection.

bp_note <- function(Reg_OECDEurope_bp) {
  paste0(
    "\"Best Practice\" predictions are derived from the RoBMA meta-regression in Table 10. ",
    "Predictions are evaluated at SC1_Cognitive = 1, LaggedDV = 0, NumberSCVars = 1, Endog_FE = 1, ",
    "Reg_OECDEurope = ", Reg_OECDEurope_bp, ", Reg_US = 0, Reg_Africa = 0, Reg_Asia = 0 ",
    "(Reg_Other implied as the reference category when all four region dummies are 0). ",
    "se(z) and publication year are not covariates in this model: RoBMA corrects for publication ",
    "bias internally rather than through a linear se(z) term, and publication year was not retained ",
    "by the BIC pre-selection in Table 9."
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
# Three sheets, matching Part 1 Table 7's structure exactly: one panel
# per best-practice scenario, plus a Notes sheet with the full covariate
# list for each.

write_xlsx(
  list(
    "BP1 - Non-OECD Europe" = panel_bp1,
    "BP2 - OECD Europe"     = panel_bp2,
    "Notes"                 = notes_sheet
  ),
  "Table11_RoBMA_BestPractice.xlsx"
)
cat("\nResults saved to Table11_RoBMA_BestPractice.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds\n", total_time["elapsed"]))
cat("\nProtocol_Table11.R complete.\n")
