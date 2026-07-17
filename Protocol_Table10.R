library(haven)
library(RoBMA)
library(openxlsx)

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record overall start time
start_time <- proc.time()

# ---- LOAD DATA ----
dat <- readRDS("SCData_processed.rds")
cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- FISHER'S Z TRANSFORMATION ----
# z = atanh(pcc); vi = 1 / (df - 1)
dat$yi <- atanh(dat$pcc)
dat$vi <- 1 / (dat$df - 1)

# ---- VARIABLE TYPES ----
# Binary (0/1) indicators: kept as integer. For a binary predictor, numeric
# and a two-level factor produce identical model matrices, so factor conversion
# is unnecessary. Factor conversion IS required for multi-level categoricals
# (e.g., CityLevel) to prevent R from treating them as continuous.
binary_mods <- c("SC1_Cognitive", "LaggedDV", "Endog_FE",
                 "Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia")
for (v in binary_mods) {
  dat[[v]] <- as.integer(dat[[v]])
}
cat(sprintf("\nBinary moderators declared as integer: %s\n",
            paste(binary_mods, collapse = ", ")))

# Continuous moderator: passed in raw form. RoBMA automatically standardizes
# continuous predictors internally (standardize_continuous_predictors = TRUE
# by default), making the prior scale-invariant and centering the intercept
# on the grand mean effect. No manual standardization is needed, and
# Protocol_Table11.R can pass raw covariate values directly to predict().
cat(sprintf("NumberSCVars passed as raw count (RoBMA standardizes internally).\n\n"))

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
# Selected via BIC drop-one in Table 9 (Protocol_Table9.R), plus the four
# region dummies added to enable best-practice predictions for two regional
# scenarios in Table 11.
#
#   SC1_Cognitive   -- retained (convergence failure when dropped = retain by default)
#   LaggedDV        -- retained (convergence failure when dropped = retain by default)
#   NumberSCVars    -- retained (delta-BIC = +1.42); continuous, standardized internally by RoBMA
#   Endog_FE        -- retained (delta-BIC = +1.05)
#   Reg_OECDEurope  -- added to support Table 11 best-practice predictions
#   Reg_US          -- added to support Table 11 best-practice predictions
#   Reg_Africa      -- added to support Table 11 best-practice predictions
#   Reg_Asia        -- added to support Table 11 best-practice predictions
#
# The four region dummies (Reg_OECDEurope, Reg_US, Reg_Africa, Reg_Asia) are
# included as a set even though none passed the BIC threshold in Table 9. This
# mirrors Part 1 (Table 7), where all four are included with Reg_Other as the
# omitted reference category. Including all four is necessary so that Table 11
# can correctly specify BP#1 (Non-OECD/Europe: all four = 0, implying Reg_Other)
# and BP#2 (OECD/Europe: Reg_OECDEurope = 1, others = 0). Without all four,
# setting Reg_OECDEurope = 0 would predict for a mix of US, Africa, Asia, and
# Other studies rather than specifically Other.

# ---- FIT MULTILEVEL RoBMA META-REGRESSION ----
# cluster = newid: effect sizes nested within studies (multilevel structure)
# measure = "ZCOR": input is Fisher's z transformed correlations
# mods: substantive moderators only (no sez -- see note above)
#
# NOTE ON COMPUTATION TIME:
# The most recent full run of this script took approximately 4 hours
# (14,400 seconds) on standard hardware.

cat("Fitting multilevel RoBMA meta-regression...\n")
cat("Expected runtime: ~4 hours based on prior run.\n\n")

time_fit <- system.time({
  fit <- RoBMA(
    yi       = yi,
    vi       = vi,
    measure  = "ZCOR",
    mods     = ~ SC1_Cognitive + LaggedDV + NumberSCVars + Endog_FE +
                 Reg_OECDEurope + Reg_US + Reg_Africa + Reg_Asia,
    cluster  = newid,
    sample   = 20000,
    burnin   = 10000,
    adapt    = 10000,
    thin     = 5,
    parallel = TRUE,
    seed     = 1,
    data     = dat
  )
})

cat(sprintf("\nModel fitting time: %.1f seconds (%.1f hours)\n",
            time_fit["elapsed"], time_fit["elapsed"] / 3600))

# ---- SAVE MODEL OBJECT ----
# Protocol_Table11.R loads this object directly to avoid refitting.
cat("Saving model object...\n")
saveRDS(fit, "RoBMA_metaregression.rds")
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

# ---- EXPORT TO EXCEL ----
# Five sheets: component inclusion, moderator inclusion, regression
# coefficients, common estimates (tau/rho), publication bias weights.

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

saveWorkbook(wb, "Table10_RoBMA_MetaRegression.xlsx", overwrite = TRUE)
cat("Results saved to Table10_RoBMA_MetaRegression.xlsx\n")

# ---- OVERALL RUN TIME ----
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f hours)\n",
            total_time["elapsed"], total_time["elapsed"] / 3600))

cat("\nProtocol_Table10.R complete. Proceed to Protocol_Table11.R for best-practice predictions.\n")
