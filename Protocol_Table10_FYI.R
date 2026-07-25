library(here)     # resolves file paths relative to the project root (.Rproj)
library(haven)
library(RoBMA)
library(openxlsx)

# Record overall start time
start_time <- proc.time()

# ============================================================
# FYI / ROBUSTNESS CHECK -- NOT THE PRIMARY PROTOCOL RESULT
# ============================================================
# Protocol_Table10.R implements the protocol's actual decision rule: the
# moderators retained by Table 9's BIC drop-one pre-selection (currently
# just Endog_FE), plus Reg_OECDEurope for Table 11's best-practice
# scenarios. That is the two-moderator model this protocol reports as its
# primary RoBMA meta-regression.
#
# This script asks a different, narrower question: how sensitive is the
# Effect-component verdict to that specific moderator choice? Rather than
# BIC-selected moderators, it uses the three moderators Table 6's CHE model
# itself classifies as "Evidence of Moderator Effect" in Panel B --
# Endogeneity: FE, Endogeneity: IV, and Lagged dependent variable (CHE
# p < 0.10) -- plus Reg_OECDEurope, again to support best-practice
# predictions in Protocol_Table11_FYI.R.
#
# Sparsity check before fitting (see CLAUDE.md for the full comparison):
# Endog_FE (18 studies) and LaggedDV (9 studies) have already been shown to
# converge cleanly in the primary Table 10 model or its five-moderator
# predecessor; Endog_IV (29 studies) is new to a RoBMA fit but is better
# populated than Reg_OECDEurope (27 studies), which itself converged
# (R-hat ~1.04) in every run so far. None of the four moderators here are
# as sparse as the region dummies (Reg_Asia: 8 studies, Reg_US: 6, Reg_
# Africa: 4) that caused severe non-convergence in an earlier specification.
#
# UPDATE after the first run (session 14): the sparsity assessment above
# was too optimistic. summary(fit, diagnostics = TRUE) on that run showed
# the intercept (R-hat 1.051, ESS 238), LaggedDV (R-hat 1.080, ESS 90 --
# the worst diagnostic seen anywhere in this project), and Reg_OECDEurope
# (ESS 339) all failing the protocol's own bar (R-hat < 1.05, ESS >= 500),
# even with autofit = TRUE already on. LaggedDV, the sparsest moderator
# here (9 studies), is the likely driver -- the same product-space-MCMC
# entanglement mechanism diagnosed for the sparse region dummies earlier
# in this project (see Protocol_Table10.R) -- and it appears to have
# dragged the intercept's own mixing down with it, which directly
# undermines the Effect-component Bayes factor (0.460 in that run) since
# the Effect BF is computed from the intercept. Per the user's decision,
# sample/burnin/adapt below have been doubled (autofit alone did not fix
# it last time) and the model is being refit from scratch before any of
# these numbers are treated as reliable.
#
# Results from this script are reported alongside the primary Table 10/11
# results as a robustness check, not as a replacement for them.

# ---- LOAD DATA ----
dat <- readRDS(here("SCData_processed.rds"))
cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- FISHER'S Z TRANSFORMATION ----
# z = atanh(pcc); vi = 1 / (df - 1)
dat$yi <- atanh(dat$pcc)
dat$vi <- 1 / (dat$df - 1)

# ---- VARIABLE TYPES ----
# Binary (0/1) indicators: declared as factors, matching Protocol_Table10.R's
# convention and for the same reason (RoBMA applies the categorical-moderator
# prior to factors and the continuous-predictor prior, on the wrong scale, to
# 0/1 integers).
binary_mods <- c("Endog_FE", "Endog_IV", "LaggedDV", "Reg_OECDEurope")
for (v in binary_mods) {
  dat[[v]] <- as.factor(dat[[v]])
}
cat(sprintf("\nBinary moderators declared as factors: %s\n",
            paste(binary_mods, collapse = ", ")))

# ---- WHY sez IS NOT INCLUDED AS A MODERATOR ----
# Same reasoning as Protocol_Table10.R: RoBMA corrects for publication bias
# internally through its model-averaging framework, so sez is not entered as
# an explicit moderator here either.

# ---- WHY ONLY ONE REGION DUMMY (not four) ----
# Same reasoning as Protocol_Table10.R: neither best-practice scenario in
# Protocol_Table11_FYI.R distinguishes among Reg_US, Reg_Africa, and
# Reg_Asia, and estimating them separately caused severe non-convergence in
# an earlier specification (see Protocol_Table10.R and CLAUDE.md for the
# full diagnosis). Those three remain folded into the reference category
# here as well.

# ---- FIT MULTILEVEL RoBMA META-REGRESSION ----
# cluster = newid: effect sizes nested within studies (multilevel structure)
# measure = "ZCOR": input is Fisher's z transformed correlations
# mods: the three CHE-significant moderators from Table 6 Panel B, plus
# Reg_OECDEurope for best-practice scenario predictions.
#
# NOTE ON COMPUTATION TIME:
# The first attempt at this fit (standard budget + autofit = TRUE) took 3.5
# hours and still failed convergence on the intercept, LaggedDV, and
# Reg_OECDEurope (see the UPDATE note above). sample/burnin/adapt are
# doubled below, matching the fix that worked for the primary Table 10
# model's region-dummy convergence problem in session 9. Expect roughly
# double the previous runtime -- perhaps 6-8 hours -- possibly less, but
# treat that as a rough expectation rather than a guarantee.

cat("Fitting multilevel RoBMA meta-regression (FYI robustness check)...\n")
cat("Expected runtime: roughly 6-8 hours (doubled MCMC budget after the first attempt failed convergence).\n\n")

time_fit <- system.time({
  fit <- RoBMA(
    yi       = yi,
    vi       = vi,
    measure  = "ZCOR",
    mods     = ~ Endog_FE + Endog_IV + LaggedDV + Reg_OECDEurope,
    cluster  = newid,
    # Factor moderators keep RoBMA's default mean-difference ("meandif")
    # contrasts, matching Protocol_Table10.R, so the model-averaged
    # intercept remains the grand-mean (average-study) effect.
    sample   = 40000,   # doubled from the protocol-standard 20000 after
    burnin   = 20000,   # the first attempt (standard budget + autofit)
    adapt    = 20000,   # failed convergence on intercept/LaggedDV/
                        # Reg_OECDEurope -- see UPDATE note above.
    thin     = 5,
    parallel = TRUE,
    seed     = 1,
    autofit  = TRUE,    # extends chains automatically if convergence
                        # criteria are not met within the (now doubled) budget
    data     = dat
  )
})

cat(sprintf("\nModel fitting time: %.1f seconds (%.1f hours)\n",
            time_fit["elapsed"], time_fit["elapsed"] / 3600))

# ---- SAVE MODEL OBJECT ----
# Protocol_Table11_FYI.R loads this object directly to avoid refitting.
cat("Saving model object...\n")
saveRDS(fit, here("RoBMA_metaregression_FYI.rds"))
cat("Model saved to RoBMA_metaregression_FYI.rds\n\n")

# ---- SUMMARY ----
# Reports model-averaged estimates for mu and each moderator coefficient,
# inclusion Bayes factors for effect, heterogeneity, and publication bias,
# AND full R-hat/ESS diagnostics (diagnostics = TRUE), printed directly
# this time so convergence can be checked from this run's own console
# output rather than requiring a separate interactive call. Bar: R-hat <
# 1.05, ESS >= 500 per parameter (same standard as Protocol_Table10.R).
s <- summary(fit, diagnostics = TRUE)
cat("Model summary (with MCMC diagnostics):\n")
print(s)

# Diagnostic: print field names of the summary object. If the Excel export
# below fails, check these names and update the field references accordingly.
cat("\nSummary object fields:", paste(names(s), collapse = ", "), "\n\n")

# ---- AVERAGE-STUDY BIAS-CORRECTED MEAN (for comparison to Table 12) ----
# Same construction as Protocol_Table10.R: under RoBMA's meandif contrasts,
# the intercept marginal mean is the model-averaged, publication-bias-
# corrected mean for an average study. Reported here for direct comparison
# to Table 12's "unconditional" RoBMA row, not as a replacement for it.
cat("Computing average-study marginal mean...\n")
mm       <- marginal_means(fit)
mu_draws <- as.numeric(mm$inference$averaged$mu_intercept$intercept)
avg_study <- data.frame(
  Quantity = "Average-study bias-corrected mean (Fisher's z) -- FYI model",
  Mean     = mean(mu_draws),
  `2.5%`   = unname(quantile(mu_draws, 0.025)),
  `97.5%`  = unname(quantile(mu_draws, 0.975)),
  check.names = FALSE
)
cat(sprintf("Average-study bias-corrected mean: %.3f [%.3f, %.3f]\n\n",
            avg_study$Mean, avg_study$`2.5%`, avg_study$`97.5%`))

# ---- EXPORT TO EXCEL ----
# Same six-sheet structure as Table10_RoBMA_MetaRegression.xlsx, saved under
# a distinct filename so it does not overwrite the primary protocol result.

# Helper: convert a summary table (with row names) to a plain data frame
to_df <- function(x, label_col) {
  df <- as.data.frame(x)
  df <- cbind(setNames(data.frame(rownames(df), stringsAsFactors = FALSE), label_col), df)
  rownames(df) <- NULL
  df
}

wb <- createWorkbook()

# Sheet 1: Component inclusion (Effect, Heterogeneity, Publication Bias)
addWorksheet(wb, "ComponentInclusion")
writeData(wb, "ComponentInclusion",
          to_df(s$inclusion_components, "Component"))

# Sheet 2: Moderator inclusion (PIPs and BFs for each moderator)
addWorksheet(wb, "ModeratorInclusion")
writeData(wb, "ModeratorInclusion",
          to_df(s$inclusion_mods, "Moderator"))

# Sheet 3: Model-averaged regression coefficients (intercept + moderators)
addWorksheet(wb, "Coefficients")
writeData(wb, "Coefficients",
          rbind(to_df(s$estimates_scale, "Parameter"),
                to_df(s$estimates_mods,  "Parameter")))

# Sheet 4: Common estimates (tau, rho)
addWorksheet(wb, "CommonEstimates")
writeData(wb, "CommonEstimates",
          to_df(s$estimates, "Parameter"))

# Sheet 5: Publication bias weights
addWorksheet(wb, "PublicationBias")
writeData(wb, "PublicationBias",
          to_df(s$estimates_bias, "Parameter"))

# Sheet 6: Average-study bias-corrected mean (FYI model)
addWorksheet(wb, "AverageStudy")
writeData(wb, "AverageStudy", avg_study)

saveWorkbook(wb, here("Table10_FYI_RoBMA.xlsx"), overwrite = TRUE)
cat("Results saved to Table10_FYI_RoBMA.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f hours)\n",
            total_time["elapsed"], total_time["elapsed"] / 3600))

cat("\nProtocol_Table10_FYI.R complete. Proceed to Protocol_Table11_FYI.R for best-practice predictions.\n")
