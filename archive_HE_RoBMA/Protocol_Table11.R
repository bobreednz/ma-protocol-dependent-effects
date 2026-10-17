library(here)   # resolves file paths relative to the project root (.Rproj)
library(RoBMA)
library(writexl)

start_time <- proc.time()

# ---- LOAD MODEL ----
# The fitted multilevel RoBMA meta-regression from Protocol_Table10.R. No refitting needed.
cat("Loading Protocol_Table10.R model...\n")
fit <- readRDS(here("RoBMA_metaregression.rds"))
cat("Model loaded.\n\n")

# ---- BEST-PRACTICE COVARIATE VALUES ----
# Following Part 1 (Table 7) and Xue, Reed, and van Aert (2025), the
# best-practice study uses fixed effects for endogeneity. sez and PubYear are
# not in the model (RoBMA corrects for publication bias internally; PubYear
# was not retained by Table 9's BIC pre-selection).
#
# As of the 2026-07-24 respecification (see Protocol_Table10.R and
# CLAUDE.md), Table 10 carries only two moderators: Endog_FE and one region
# dummy, Reg_OECDEurope. SC1_Cognitive, LaggedDV, and NumberSCVars were
# dropped -- a 2026-07-23 re-run of Table 9's BIC pre-selection, independently
# reproduced 2026-07-24, no longer retains them (the four-moderator set used
# in an earlier version of this script came from a stale, pre-2026-07-15
# Table 9 run; see CLAUDE.md). Reg_US, Reg_Africa, and Reg_Asia remain folded
# into the reference category, for the reasons given in Protocol_Table10.R:
# neither best-practice scenario below distinguishes among them, and
# estimating them separately caused severe MCMC non-convergence.
#   BP#1: Non-OECD/Europe (Reg_OECDEurope = 0, implying the combined reference
#         category of Reg_US, Reg_Africa, Reg_Asia, and Reg_Other)
#   BP#2: OECD/Europe     (Reg_OECDEurope = 1)
#
# The binary moderator is encoded as a factor with levels c(0, 1) to match the
# factor coding used when the model was fit in Protocol_Table10.R, so that
# predict() applies the same contrasts.

bp1 <- data.frame(
  Endog_FE       = factor(1, levels = c(0, 1)),
  Reg_OECDEurope = factor(0, levels = c(0, 1))   # Non-OECD/Europe
)

bp2 <- data.frame(
  Endog_FE       = factor(1, levels = c(0, 1)),
  Reg_OECDEurope = factor(1, levels = c(0, 1))   # OECD/Europe
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
# RoBMA's meta-regression (Table 10) uses the moderators retained after BIC
# pre-selection (Table 9) plus one region dummy, Reg_OECDEurope (Reg_US,
# Reg_Africa, and Reg_Asia are folded into the reference category -- see
# Protocol_Table10.R). Unlike Part 1 Table 7, se(z) and publication year are
# not covariates here: RoBMA corrects for publication bias internally rather
# than through a linear se(z) term, and publication year was not among the
# moderators retained by Table 9's BIC pre-selection. SC1_Cognitive, LaggedDV,
# and NumberSCVars are also not covariates here as of the 2026-07-24
# respecification -- the current, reproducible Table 9 BIC run does not
# retain them (see Protocol_Table10.R and CLAUDE.md).

bp_note <- function(Reg_OECDEurope_bp) {
  paste0(
    "\"Best Practice\" predictions are derived from the RoBMA meta-regression in Table 10. ",
    "Predictions are evaluated at Endog_FE = 1, ",
    "Reg_OECDEurope = ", Reg_OECDEurope_bp, " ",
    "(the reference category -- combining Reg_US, Reg_Africa, Reg_Asia, and Reg_Other -- is implied ",
    "when Reg_OECDEurope = 0). se(z) and publication year are not covariates in this model: RoBMA ",
    "corrects for publication bias internally rather than through a linear se(z) term, and publication ",
    "year was not retained by the BIC pre-selection in Table 9. No other moderator is retained by ",
    "that pre-selection, so the model conditions on Endog_FE and Reg_OECDEurope only (see the ",
    "Decision point following Table 10)."
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
  here("Table11_RoBMA_BestPractice.xlsx")
)
cat("\nResults saved to Table11_RoBMA_BestPractice.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds\n", total_time["elapsed"]))
cat("\nProtocol_Table11.R complete.\n")
