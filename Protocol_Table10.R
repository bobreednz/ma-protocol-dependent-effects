library(here)     # resolves file paths relative to the project root (.Rproj)
library(haven)
library(RoBMA)
library(openxlsx)

# Record overall start time
start_time <- proc.time()

# ---- LOAD DATA ----
dat <- readRDS(here("SCData_processed.rds"))
cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- FISHER'S Z TRANSFORMATION ----
# z = atanh(pcc); vi = 1 / (df - 1)
dat$yi <- atanh(dat$pcc)
dat$vi <- 1 / (dat$df - 1)

# ---- VARIABLE TYPES ----
# Binary (0/1) indicators: declared as factors. Unlike OLS or metafor, where a
# 0/1 numeric and a two-level factor give identical fits, RoBMA treats the two
# types differently. A numeric or integer predictor is treated as continuous and
# internally standardized (standardize_predictors = TRUE by default), which
# applies the continuous-predictor prior on the wrong (per-SD) scale for a binary
# indicator. A factor is instead given the categorical-moderator prior with mean-
# difference ("meandif") contrasts, which is RoBMA's default. We keep that default
# rather than switching to treatment (dummy) contrasts; the Section 9 discussion
# explains why. The binaries are therefore declared as factors here.
binary_mods <- c("Endog_FE", "Reg_OECDEurope")
for (v in binary_mods) {
  dat[[v]] <- as.factor(dat[[v]])
}
cat(sprintf("\nBinary moderators declared as factors: %s\n",
            paste(binary_mods, collapse = ", ")))

# ---- WHY sez IS NOT INCLUDED AS A MODERATOR ----
# In FAT-PET (Part 1 Table 5/6), sez must be included explicitly because the
# only way to correct for publication bias is to set SE = 0 in the regression.
# RoBMA is different: it corrects for publication bias internally, through a
# model-averaging framework that includes 6 selection function models plus
# PET and PEESE sub-models. This correction operates regardless of what
# substantive moderators are included. Adding sez as an explicit moderator
# would double-count the precision effect. The moderators below are therefore
# purely substantive study characteristics.

# ---- MODERATORS ----
# Selected via BIC drop-one in Table 9 (Protocol_Table9.R), plus one region
# dummy added to enable best-practice predictions for two regional scenarios
# in Table 11.
#
#   Endog_FE        -- retained (delta-BIC = +1.69)
#   Reg_OECDEurope  -- added to support Table 11 best-practice predictions
#
# ---- WHY THIS DIFFERS FROM AN EARLIER MODERATOR SET ----
# An earlier run of Protocol_Table9.R (2026-07-03) retained four moderators
# instead of one -- SC1_Cognitive and LaggedDV via convergence-failure-
# default, NumberSCVars (delta-BIC +1.42), and Endog_FE (delta-BIC +1.05) --
# and an earlier version of this script was fit with that five-covariate
# set. A 2026-07-23 re-run of Table 9, independently reproduced on
# 2026-07-24 (delta-BIC values matching to 2-3 decimals), consistently
# retains only Endog_FE, with clean convergence on all 18 drop-one fits.
# The discrepancy was traced to a local metafor package update on
# 2026-07-15 -- squarely between the two Table 9 runs that disagree -- not
# to any data or code change: the 2026-07-03 and 2026-07-23 versions of
# Protocol_Table9.R are identical apart from here()-path mechanics. The
# pre-update four-moderator set is therefore treated as stale, and this
# script has been refit using only what the current, reproducible Table 9
# run supports. See CLAUDE.md for the full investigation.
#
# ---- WHY ONLY ONE REGION DUMMY (not four) ----
# Part 1 (Table 7) includes all four region dummies (Reg_OECDEurope, Reg_US,
# Reg_Africa, Reg_Asia, with Reg_Other as the omitted reference), because CHE
# is a frequentist GLS model with no MCMC estimation cost to including
# categories it doesn't strictly need. RoBMA's meta-regression is different:
# each categorical moderator adds a discrete inclusion/exclusion dimension to
# a single joint (product-space) MCMC chain, and Table 11's two best-practice
# scenarios never actually need to distinguish among Reg_US, Reg_Africa, and
# Reg_Asia -- BP#1 (non-OECD/Europe) sets all of them to 0, and BP#2 (OECD/
# Europe) only ever sets Reg_OECDEurope = 1, holding the rest at 0. Estimating
# three additional categories the predictions never separate is unnecessary
# complexity, and here it is also actively harmful: on 2026-07-23, fitting all
# four region dummies produced severe non-convergence specifically on the
# three sparsest ones (Reg_Asia: 8 of 83 studies carry it, R-hat 1.142, ESS
# 126; Reg_US: 6 studies, ESS 98; Reg_Africa: 4 studies) that persisted across
# a doubled sampling budget and did not resolve with autofit = TRUE. A vif()
# check ruled out ordinary collinearity (all VIFs < 1.3) but showed high
# posterior correlation between each sparse dummy's inclusion indicator and
# its own coefficient (e.g. 0.765 for Reg_Asia), and among the sparse dummies
# themselves -- consistent with known difficulty in product-space Bayesian
# model averaging when a categorical level is supported by very few clusters.
#
# This is not a fix specific to this dataset: as a general rule for this
# protocol, a categorical moderator should not be split more finely in a
# RoBMA meta-regression than the analysis actually requires, and any level
# supported by only a handful of clusters is a convergence risk worth
# checking for (via summary(fit, diagnostics = TRUE) and vif()) before
# committing to a specification, regardless of what a corresponding
# frequentist model can absorb for free. Here, since neither best-practice
# scenario in Table 11 distinguishes Reg_US, Reg_Africa, or Reg_Asia from
# each other, those three are folded into the reference category (along with
# Reg_Other) and only Reg_OECDEurope is retained as a moderator. See
# CLAUDE.md for the full non-convergence diagnosis this decision is based on.

# ---- FIT MULTILEVEL RoBMA META-REGRESSION ----
# cluster = newid: effect sizes nested within studies (multilevel structure)
# measure = "ZCOR": input is Fisher's z transformed correlations
# mods: substantive moderators only (no sez -- see note above)
#
# NOTE ON COMPUTATION TIME:
# Prior full runs with five moderators and the protocol-standard sampling
# budget took approximately 3-4 hours on standard hardware. With only two
# moderators (Endog_FE, Reg_OECDEurope) in this respecification, this run
# should be somewhat faster, but treat that as a rough expectation rather
# than a guarantee.

cat("Fitting multilevel RoBMA meta-regression...\n")
cat("Expected runtime: roughly 3-4 hours, possibly less.\n\n")

time_fit <- system.time({
  fit <- RoBMA(
    yi       = yi,
    vi       = vi,
    measure  = "ZCOR",
    mods     = ~ Endog_FE + Reg_OECDEurope,
    cluster  = newid,
    # Factor moderators keep RoBMA's default mean-difference ("meandif")
    # contrasts, so the model-averaged intercept is the grand-mean (average-study)
    # effect. This is used as the "unconditional" RoBMA corrected mean in Table 12
    # (see the marginal_means() step below and the Section 9 discussion).
    sample   = 20000,   # protocol-standard settings (matches Protocol_Figure3.R
    burnin   = 10000,   # and the "Standard MCMC settings" note in the RoBMA
    adapt    = 10000,   # Technical Notes), not the doubled budget tried during
    thin     = 5,       # troubleshooting. The doubled budget did not resolve
    parallel = TRUE,    # non-convergence while all four region dummies were in
    seed     = 1,       # the model, and much of Reg_OECDEurope's own borderline
    autofit  = TRUE,    # ESS (354, just under the 500 threshold) may have come
                        # from being entangled with the sparser dummies now
                        # removed (posterior correlations of 0.15-0.23 with
                        # Reg_Asia/Reg_US -- see CLAUDE.md). Reverting to the
                        # standard budget is a genuine test of whether the
                        # specification fix, not extra compute, was the actual
                        # solution; keeping a doubled budget as a standing
                        # exception would also make the "standard MCMC
                        # settings" claim in the protocol text inaccurate.
    data     = dat
  )
})

cat(sprintf("\nModel fitting time: %.1f seconds (%.1f hours)\n",
            time_fit["elapsed"], time_fit["elapsed"] / 3600))

# ---- SAVE MODEL OBJECT ----
# Protocol_Table11.R loads this object directly to avoid refitting.
cat("Saving model object...\n")
saveRDS(fit, here("RoBMA_metaregression.rds"))
cat("Model saved to RoBMA_metaregression.rds\n\n")

# ---- SUMMARY ----
# Reports model-averaged estimates for mu and each moderator coefficient,
# plus inclusion Bayes factors for effect, heterogeneity, and publication bias.
s <- summary(fit, include_mcmc_diagnostics = FALSE)
cat("Model summary:\n")
print(s)

# Diagnostic: print field names of the summary object. If the Excel export
# below fails, check these names and update the field references accordingly.
cat("\nSummary object fields:", paste(names(s), collapse = ", "), "\n\n")

# ---- AVERAGE-STUDY BIAS-CORRECTED MEAN (for Table 12) ----
# marginal_means() returns model-averaged marginal means. Under RoBMA's default
# mean-difference ("meandif") contrasts, the "intercept" marginal mean is the
# overall grand-mean effect: the model-averaged (publication-bias-corrected) mean
# for an average study, with the continuous moderator at its mean and the factor
# moderators averaged symmetrically across their levels. This is the RoBMA
# analogue of the FAT-PET Panel B intercept (Table 5), which is evaluated at the
# moderator means, and it is carried into Table 12 as the "unconditional" RoBMA
# corrected mean. (That this value sits well below the raw mean of ~0.175, close
# to the intercept-only bias-corrected mean in Table 8, confirms the marginal
# mean is bias-corrected.) The credible interval is computed directly from the
# model-averaged posterior draws for the intercept marginal mean.
cat("Computing average-study marginal mean...\n")
mm       <- marginal_means(fit)
mu_draws <- as.numeric(mm$inference$averaged$mu_intercept$intercept)
avg_study <- data.frame(
  Quantity = "Average-study bias-corrected mean (Fisher's z)",
  Mean     = mean(mu_draws),
  `2.5%`   = unname(quantile(mu_draws, 0.025)),
  `97.5%`  = unname(quantile(mu_draws, 0.975)),
  check.names = FALSE
)
cat(sprintf("Average-study bias-corrected mean: %.3f [%.3f, %.3f]\n\n",
            avg_study$Mean, avg_study$`2.5%`, avg_study$`97.5%`))

# ---- EXPORT TO EXCEL ----
# Six sheets: component inclusion, moderator inclusion, regression coefficients,
# common estimates (tau/rho), publication bias weights, and the average-study
# bias-corrected mean (for Table 12).

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
# s$estimates_scale contains the intercept; s$estimates_mods contains the
# moderator coefficients. Combine into one sheet.
addWorksheet(wb, "Coefficients")
writeData(wb, "Coefficients",
          rbind(to_df(s$estimates_scale, "Parameter"),
                to_df(s$estimates_mods,  "Parameter")))

# Sheet 4: Common estimates (tau, rho)
# s$estimates contains tau and rho (not the regression coefficients)
addWorksheet(wb, "CommonEstimates")
writeData(wb, "CommonEstimates",
          to_df(s$estimates, "Parameter"))

# Sheet 5: Publication bias weights
addWorksheet(wb, "PublicationBias")
writeData(wb, "PublicationBias",
          to_df(s$estimates_bias, "Parameter"))

# Sheet 6: Average-study bias-corrected mean (grand-mean marginal mean), the
# "unconditional" RoBMA corrected mean carried into Table 12.
addWorksheet(wb, "AverageStudy")
writeData(wb, "AverageStudy", avg_study)

saveWorkbook(wb, here("Table10_RoBMA_MetaRegression.xlsx"), overwrite = TRUE)
cat("Results saved to Table10_RoBMA_MetaRegression.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f hours)\n",
            total_time["elapsed"], total_time["elapsed"] / 3600))

cat("\nProtocol_Table10.R complete. Proceed to Protocol_Table11.R for best-practice predictions.\n")
