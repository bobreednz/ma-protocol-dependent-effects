library(here)          # resolves file paths relative to the project root (.Rproj)
library(metafor)      # rma.mv() for CHE meta-regression
library(openxlsx)     # save xlsx

# Record overall start time
start_time <- proc.time()

# ---- LOAD DATA ----
dat <- readRDS(here("SCData_processed.rds"))

# Compute Fisher's z and its standard error from the raw partial correlations.
# z = atanh(pcc); sez = sqrt(1 / (df - 1)), the SE of Fisher's z for a PCC.
dat$z   <- atanh(dat$pcc)
dat$sez <- sqrt(1 / (dat$df - 1))

cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- DEFINE THE 18 MODERATORS ----
# sez is the publication bias proxy; it is NOT a candidate for dropping --
# it is always retained regardless of BIC. We evaluate the 18 study
# characteristics only.
moderators <- c(
  "SC1_Cognitive",  # social capital type: cognitive
  "SC1_Structural", # social capital type: structural
  "DV_GrowthRate",  # outcome: GDP growth rate
  "PubYear",        # publication year (continuous)
  "Published",      # 1 = published in peer-reviewed journal
  "LaggedDV",       # 1 = lagged dependent variable included
  "LaggedSC",       # 1 = lagged social capital used
  "NumberSCVars",   # number of social capital variables
  "Endog_IV",       # 1 = instrumental variable used
  "Endog_FE",       # 1 = fixed effects used
  "CityLevel",      # 1 = city-level analysis
  "RegionLevel",    # 1 = region-level analysis
  "CountryLevel",   # 1 = country-level analysis
  "PanelData",      # 1 = panel data
  "Reg_OECDEurope", # 1 = OECD / European sample
  "Reg_US",         # 1 = US sample
  "Reg_Africa",     # 1 = African sample
  "Reg_Asia"        # 1 = Asian sample
)

# All 19 regressors (sez always included)
all_regressors <- c("sez", moderators)

# ---- BUILD IMPUTED COVARIANCE MATRIX ----
# Same rho = 0.5 assumption used throughout the protocol.
# metafor::vcalc() replaced the deprecated clubSandwich::impute_covariance_matrix().
V_mat <- metafor::vcalc(
  vi      = dat$sez^2,
  cluster = dat$newid,
  rho     = 0.5
)

# ---- FIT THE FULL MODEL WITH ML ----
# IMPORTANT: We use method = "ML" (not "REML") because BIC comparisons are
# only valid when models differ in their fixed effects (moderators). REML
# log-likelihoods are not comparable across models with different fixed
# effects structures, so ML is required here.
cat("\nFitting full CHE model with ML (19 regressors)...\n")

full_formula <- as.formula(
  paste("~", paste(all_regressors, collapse = " + "))
)

time_full <- system.time({
  fit_full <- rma.mv(
    yi     = z,
    V      = V_mat,
    mods   = full_formula,
    random = ~ 1 | newid/obs,
    data   = dat,
    method = "ML"   # ML required for valid BIC comparisons
  )
})

bic_full <- BIC(fit_full)
cat(sprintf("Full model BIC: %.4f  (time: %.1f sec)\n", bic_full, time_full["elapsed"]))

# ---- DROP-ONE BIC COMPARISONS ----
# For each of the 18 moderators, refit the model without that variable and
# compute delta-BIC = BIC(dropped) - BIC(full).
# A positive delta-BIC means BIC increased when the variable was dropped,
# i.e., the variable improves model fit -- RETAIN it.
# A negative or zero delta-BIC means dropping the variable does not hurt
# fit -- EXCLUDE it from the RoBMA meta-regression.

cat("\nRunning drop-one BIC comparisons...\n")

results <- data.frame(
  Variable   = character(length(moderators)),
  BIC_full   = numeric(length(moderators)),
  BIC_dropped = numeric(length(moderators)),
  Delta_BIC  = numeric(length(moderators)),
  Retain     = character(length(moderators)),
  stringsAsFactors = FALSE
)

for (i in seq_along(moderators)) {
  var_drop <- moderators[i]

  # Build formula omitting this one variable
  remaining  <- setdiff(all_regressors, var_drop)
  drop_formula <- as.formula(
    paste("~", paste(remaining, collapse = " + "))
  )

  # Fit the reduced model
  t_start <- proc.time()
  fit_drop <- tryCatch(
    rma.mv(
      yi     = z,
      V      = V_mat,
      mods   = drop_formula,
      random = ~ 1 | newid/obs,
      data   = dat,
      method = "ML"
    ),
    error = function(e) {
      cat(sprintf("  WARNING: Model without %s failed: %s\n", var_drop, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit_drop)) {
    bic_drop  <- NA
    delta_bic <- NA
    retain    <- "ERROR"
  } else {
    bic_drop  <- BIC(fit_drop)
    delta_bic <- bic_drop - bic_full
    # Retain if dropping the variable increases BIC (delta > 0)
    retain    <- ifelse(delta_bic > 0, "Yes", "No")
  }

  results$Variable[i]    <- var_drop
  results$BIC_full[i]    <- bic_full
  results$BIC_dropped[i] <- bic_drop
  results$Delta_BIC[i]   <- delta_bic
  results$Retain[i]      <- retain

  t_elapsed <- (proc.time() - t_start)["elapsed"]
  cat(sprintf("  %-20s  delta-BIC = %+.2f  Retain: %s  (%.1f sec)\n",
              var_drop, delta_bic, retain, t_elapsed))
}

# ---- REPORT RETAINED MODERATORS ----
retained <- results$Variable[results$Retain == "Yes"]
cat(sprintf("\n%d of 18 moderators retained (delta-BIC > 0):\n", length(retained)))
for (v in retained) cat(sprintf("  %s\n", v))

cat("\nModerators retained for the RoBMA meta-regression (Protocol_Table10.R):\n")
cat(paste(retained, collapse = ", "), "\n")
cat("(sez is not carried forward -- RoBMA corrects for publication bias internally.)\n")

# ---- SAVE RESULTS TABLE AS XLSX ----
# Nicer display labels matching Table 6 conventions
var_labels <- c(
  SC1_Cognitive  = "SC: Cognitive",
  SC1_Structural = "SC: Structural",
  DV_GrowthRate  = "Outcome: Growth rate",
  PubYear        = "Publication year",
  Published      = "Published",
  LaggedDV       = "Lagged dependent variable",
  LaggedSC       = "Lagged social capital",
  NumberSCVars   = "Number of SC variables",
  Endog_IV       = "Endogeneity: IV",
  Endog_FE       = "Endogeneity: FE",
  CityLevel      = "City-level",
  RegionLevel    = "Region-level",
  CountryLevel   = "Country-level",
  PanelData      = "Panel data",
  Reg_OECDEurope = "Region: OECD/Europe",
  Reg_US         = "Region: US",
  Reg_Africa     = "Region: Africa",
  Reg_Asia       = "Region: Asia"
)

display_table <- data.frame(
  Variable  = var_labels[results$Variable],
  Delta_BIC = round(results$Delta_BIC, 3),
  Retain    = results$Retain,
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "BIC_Selection")
writeData(wb, "BIC_Selection", display_table)

# Add a note row explaining the decision rule
note <- data.frame(
  Variable  = "Note: delta-BIC = BIC(model without variable) - BIC(full model). Retain if delta-BIC > 0.",
  Delta_BIC = NA,
  Retain    = NA
)
writeData(wb, "BIC_Selection", note, startRow = nrow(display_table) + 3, colNames = FALSE)

saveWorkbook(wb, here("Table9_BIC_ModeratorSelection.xlsx"), overwrite = TRUE)
cat("Results saved to Table9_BIC_ModeratorSelection.xlsx\n")

cat("\nProtocol_Table9.R complete.\n")
cat("Proceed to Protocol_Table10.R using the retained moderators listed above.\n")

# Overall run time
total_time <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes)\n",
            total_time["elapsed"], total_time["elapsed"] / 60))
