setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("haven", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(haven)      # read/write Stata .dta files
library(tidyverse)  # data manipulation
library(writexl)    # export to Excel

# ============================================================
# Load data
# ============================================================

DT <- as.data.frame(read_dta("SCData (20240508).dta"))

cat("==========================================================\n")
cat("DATA CLEANING REPORT\n")
cat("Dataset: SCData (20240508).dta\n")
cat(sprintf("Observations: %d   |   Studies: %d\n",
            nrow(DT), length(unique(DT$newid))))
cat("==========================================================\n\n")

# ============================================================
# CHECK 1: Missing values
#
# Meta-analysis datasets should have no missing values in any
# variable used for estimation. We check all columns.
# ============================================================

cat("--- CHECK 1: Missing values ---\n")

missing_counts <- colSums(is.na(DT))
missing_vars   <- missing_counts[missing_counts > 0]

if (length(missing_vars) == 0) {
  cat("PASS: No missing values in any variable.\n\n")
} else {
  cat("FAIL: Missing values found in the following variables:\n")
  print(missing_vars)
  cat("\n")
}

# ============================================================
# CHECK 2: Duplicate observations
#
# We check for two types of duplicates:
#   (a) Duplicate obs identifiers (the obs variable should be
#       a unique row identifier).
#   (b) Exact content duplicates: rows with the same study
#       (newid), effect size (z), standard error (sez), and
#       t-statistic (tstat). These represent estimates that
#       appear more than once in the dataset. They may be
#       legitimate (e.g., the same estimate reported under
#       multiple specifications) or data entry errors.
# ============================================================

cat("--- CHECK 2: Duplicate observations ---\n")

# (a) Duplicate obs identifiers
dup_obs <- sum(duplicated(DT$obs))
if (dup_obs == 0) {
  cat("PASS: All obs identifiers are unique.\n")
} else {
  cat(sprintf("FAIL: %d duplicate obs identifiers found.\n", dup_obs))
}

# (b) True duplicates: rows identical on ALL variables except obs.
# These represent the same estimate entered twice and one copy
# should be removed before analysis.
all_vars_except_obs <- setdiff(names(DT), "obs")

DT <- DT %>%
  group_by(across(all_of(all_vars_except_obs))) %>%
  mutate(flag_true_duplicate = as.integer(n() > 1)) %>%
  ungroup()

n_true_dup <- sum(DT$flag_true_duplicate)

if (n_true_dup == 0) {
  cat("PASS: No true duplicates (identical on all variables except obs).\n")
} else {
  cat(sprintf("FAIL: %d observations (%d pairs) are identical on all variables except obs.\n",
              n_true_dup, as.integer(n_true_dup / 2)))
  cat("      One copy of each pair should be removed before analysis.\n")
  true_dup_summary <- DT %>%
    filter(flag_true_duplicate == 1) %>%
    group_by(newid) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(newid)
  cat("      Studies affected:\n")
  print(as.data.frame(true_dup_summary))
  cat("\n")
}

# (c) Coding conflicts: rows that share newid, z, sez, and tstat
# but differ on at least one other variable. These represent the
# same statistical estimate coded with conflicting study
# characteristics -- e.g., the same regression result assigned
# to different social capital subcategories, or with inconsistent
# endogeneity treatment coding. These are NOT duplicates to be
# dropped; they indicate a data coding decision that needs to
# be resolved by the researcher.

DT <- DT %>%
  group_by(newid, z, sez, tstat) %>%
  mutate(flag_coding_conflict = as.integer(n() > 1 & flag_true_duplicate == 0)) %>%
  ungroup()

n_conflict <- sum(DT$flag_coding_conflict)

if (n_conflict == 0) {
  cat("PASS: No coding conflicts found.\n\n")
} else {
  cat(sprintf("NOTE: %d observations share newid/z/sez/tstat but differ on other variables.\n",
              n_conflict))
  cat("      These are NOT duplicates to drop. They represent the same estimate\n")
  cat("      coded with conflicting study characteristics and require researcher review.\n")
  cat("      Check the 'Coding conflicts' sheet in DataCleaning_Report.xlsx\n")
  cat("      for the specific variables that differ.\n\n")
}

# For backward compatibility, flag_duplicate marks either type
DT$flag_duplicate <- as.integer(
  DT$flag_true_duplicate == 1 | DT$flag_coding_conflict == 1
)

# ============================================================
# CHECK 3: Valid ranges for key variables
#
# We verify that all key variables fall within expected ranges.
# For binary (dummy) variables, values should be 0 or 1.
# For continuous variables, we check for plausible bounds.
# ============================================================

cat("--- CHECK 3: Valid ranges ---\n")

# Binary variables: must be 0 or 1
binary_vars <- c(
  "Published", "PanelData",
  "CityLevel", "RegionLevel", "CountryLevel", "OtherLevel",
  "LaggedSC", "LaggedDV", "Endog_IV", "Endog_FE", "NoEndogeneity",
  "DV_GrowthRate", "DV_GDPLevel",
  "SC1_Cognitive", "SC1_Structural", "SC1_Other",
  "Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia", "Reg_Other"
)

all_binary_ok <- TRUE
for (v in binary_vars) {
  bad <- DT[[v]][!DT[[v]] %in% c(0, 1)]
  if (length(bad) > 0) {
    cat(sprintf("FAIL: %s contains non-binary values: %s\n", v, paste(bad, collapse = ", ")))
    all_binary_ok <- FALSE
  }
}
if (all_binary_ok) cat("PASS: All binary variables contain only 0/1 values.\n")

# sez must be strictly positive
if (all(DT$sez > 0)) {
  cat("PASS: All sez values are strictly positive.\n")
} else {
  cat(sprintf("FAIL: %d observations have sez <= 0.\n", sum(DT$sez <= 0)))
}

# df must be strictly positive
if (all(DT$df > 0)) {
  cat("PASS: All df values are strictly positive.\n")
} else {
  cat(sprintf("FAIL: %d observations have df <= 0.\n", sum(DT$df <= 0)))
}

# z must be finite (no Inf or NaN)
if (all(is.finite(DT$z))) {
  cat("PASS: All z values are finite.\n")
} else {
  cat(sprintf("FAIL: %d observations have non-finite z values.\n", sum(!is.finite(DT$z))))
}

# PubYear: reasonable range for published social capital research
if (all(DT$PubYear >= 1990 & DT$PubYear <= 2030)) {
  cat("PASS: All PubYear values are in the range [1990, 2030].\n")
} else {
  cat(sprintf("NOTE: %d observations have PubYear outside [1990, 2030].\n",
              sum(DT$PubYear < 1990 | DT$PubYear > 2030)))
}

# NumberSCVars: must be at least 1
if (all(DT$NumberSCVars >= 1)) {
  cat("PASS: All NumberSCVars values are >= 1.\n\n")
} else {
  cat(sprintf("FAIL: %d observations have NumberSCVars < 1.\n\n",
              sum(DT$NumberSCVars < 1)))
}

# ============================================================
# CHECK 4: Mutual exclusivity of dummy variable sets
#
# Several groups of dummy variables are mutually exclusive and
# exhaustive: each observation should belong to exactly one
# category. The row sum of each group should equal 1.
# ============================================================

cat("--- CHECK 4: Mutual exclusivity of dummy variable sets ---\n")

dummy_sets <- list(
  "Social capital type (SC1)"   = c("SC1_Cognitive", "SC1_Structural", "SC1_Other"),
  "Dependent variable type (DV)"= c("DV_GrowthRate", "DV_GDPLevel"),
  "Geographic level"            = c("CityLevel", "RegionLevel", "CountryLevel", "OtherLevel"),
  "Region"                      = c("Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia", "Reg_Other"),
  "Endogeneity treatment"       = c("Endog_IV", "Endog_FE", "NoEndogeneity")
)

all_exclusive_ok <- TRUE
for (label in names(dummy_sets)) {
  vars    <- dummy_sets[[label]]
  row_sum <- rowSums(DT[, vars])
  bad     <- sum(row_sum != 1)
  if (bad == 0) {
    cat(sprintf("PASS: %s sums to 1 for all observations.\n", label))
  } else {
    cat(sprintf("FAIL: %s does not sum to 1 for %d observations.\n", label, bad))
    all_exclusive_ok <- FALSE
  }
}
cat("\n")

# ============================================================
# CHECK 5: Internal consistency of derived variables
#
# The dataset contains both z (Fisher's z) and pcc (partial
# correlation coefficient), related by z = atanh(pcc). We
# verify that this relationship holds to floating point
# precision. We also verify that tstat = pcc / sepcc.
# ============================================================

cat("--- CHECK 5: Internal consistency of derived variables ---\n")

# z vs atanh(pcc)
z_diff <- abs(DT$z - atanh(DT$pcc))
if (all(z_diff < 1e-6)) {
  cat("PASS: z = atanh(pcc) holds to within 1e-6 for all observations.\n")
} else {
  cat(sprintf("FAIL: z != atanh(pcc) for %d observations (max diff = %.2e).\n",
              sum(z_diff >= 1e-6), max(z_diff)))
}

# tstat vs pcc/sepcc
tstat_diff <- abs(DT$tstat - DT$pcc / DT$sepcc)
if (all(tstat_diff < 1e-4)) {
  cat("PASS: tstat = pcc/sepcc holds to within 1e-4 for all observations.\n\n")
} else {
  cat(sprintf("FAIL: tstat != pcc/sepcc for %d observations (max diff = %.2e).\n\n",
              sum(tstat_diff >= 1e-4), max(tstat_diff)))
}

# ============================================================
# CHECK 6: Standard error consistency with degrees of freedom
#
# For a partial correlation coefficient with df degrees of
# freedom, the standard error of Fisher's z is:
#   sez = 1 / sqrt(df - 1)
#
# We compute the "SE deficit" = (1/sqrt(df-1)) - sez, i.e.,
# how much smaller the reported sez is than the theoretical
# minimum for the given df. A positive deficit means the SE
# is smaller than the formula allows.
#
# Small deficits (< 0.01) typically reflect rounding: many
# researchers assign a single representative sez to all
# estimates from the same regression table, even when df
# varies slightly across those estimates. This is common
# practice and not a data quality concern.
#
# Larger deficits (>= 0.01) indicate a more substantial
# discrepancy. The most likely explanation is that sez was
# computed using a formula corresponding to a larger sample
# than df reflects (e.g., using bivariate rather than partial
# correlation degrees of freedom, or using total N rather than
# N-k-1). These observations warrant closer inspection.
# ============================================================

cat("--- CHECK 6: Standard error consistency with degrees of freedom ---\n")
cat("    Formula: theoretical minimum sez = 1/sqrt(df - 1)\n")
cat("    SE deficit = (1/sqrt(df-1)) - sez\n\n")

DT$sez_theoretical <- 1 / sqrt(DT$df - 1)
DT$se_deficit      <- DT$sez_theoretical - DT$sez

# Flag observations by deficit size
# Minor: likely rounding artifact (0.001 < deficit <= 0.010)
# Major: meaningful discrepancy, warrants review (deficit > 0.010)
DT$flag_se_minor <- as.integer(DT$se_deficit >  0.001 & DT$se_deficit <= 0.010)
DT$flag_se_major <- as.integer(DT$se_deficit >  0.010)

n_se_ok    <- sum(DT$se_deficit <= 0.001)
n_se_minor <- sum(DT$flag_se_minor)
n_se_major <- sum(DT$flag_se_major)

cat(sprintf("  No issue (deficit <= 0.001):    %d observations\n", n_se_ok))
cat(sprintf("  Minor deficit (0.001 to 0.010): %d observations across %d studies\n",
            n_se_minor, length(unique(DT$newid[DT$flag_se_minor == 1]))))
cat(sprintf("  Major deficit (> 0.010):        %d observations across %d studies\n",
            n_se_major, length(unique(DT$newid[DT$flag_se_major == 1]))))

if (n_se_major > 0) {
  cat("\n  Observations with major SE deficit (> 0.010):\n")
  major_df <- DT[DT$flag_se_major == 1,
                  c("newid", "obs", "z", "sez", "df", "sez_theoretical", "se_deficit")]
  major_df$se_deficit      <- round(major_df$se_deficit,      4)
  major_df$sez_theoretical <- round(major_df$sez_theoretical, 4)
  print(major_df, row.names = FALSE)
}
cat("\n")

# ============================================================
# CHECK 7: Outlier effect sizes
#
# We flag observations with unusually large effect sizes.
# In the social capital / economic growth literature, partial
# correlations above 0.80 (|z| > 1.0) are exceptionally
# large. We also flag based on a statistical criterion:
# estimates more than 3 standard deviations from the mean z.
#
# Note: flagging does not imply deletion. Large effects may
# reflect genuine heterogeneity or influential studies, and
# should be investigated before any decision to exclude.
# ============================================================

cat("--- CHECK 7: Outlier effect sizes ---\n")

# Threshold 1: |z| > 1.0 (|pcc| > 0.762)
DT$flag_outlier_z1 <- as.integer(abs(DT$z) > 1.0)

# Threshold 2: 3 standard deviations from the mean of z
z_mean <- mean(DT$z)
z_sd   <- sd(DT$z)
DT$flag_outlier_z2 <- as.integer(abs(DT$z - z_mean) > 3 * z_sd)

cat(sprintf("  Mean z = %.3f,  SD z = %.3f\n", z_mean, z_sd))
cat(sprintf("  3 SD threshold: z outside [%.3f, %.3f]\n",
            z_mean - 3 * z_sd, z_mean + 3 * z_sd))

n_out1 <- sum(DT$flag_outlier_z1)
n_out2 <- sum(DT$flag_outlier_z2)

cat(sprintf("  |z| > 1.0 (|pcc| > 0.762):  %d observations\n", n_out1))
cat(sprintf("  |z| > mean +/- 3 SD:         %d observations\n", n_out2))

if (n_out1 > 0) {
  cat("\n  Observations with |z| > 1.0:\n")
  out_df <- DT[DT$flag_outlier_z1 == 1,
                c("newid", "obs", "z", "pcc", "sez", "df", "tstat")]
  out_df[, c("z", "pcc", "sez", "tstat")] <-
    round(out_df[, c("z", "pcc", "sez", "tstat")], 4)
  print(out_df, row.names = FALSE)
}
cat("\n")

# ============================================================
# CHECK 8: Very small degrees of freedom
#
# Estimates based on very few observations are inherently
# imprecise and may be influential. We flag observations with
# df < 10, which corresponds to regressions with at most
# roughly 11-15 original observations depending on the number
# of regressors.
# ============================================================

cat("--- CHECK 8: Very small degrees of freedom ---\n")

DT$flag_small_df <- as.integer(DT$df < 10)
n_small_df       <- sum(DT$flag_small_df)

if (n_small_df == 0) {
  cat("PASS: No observations with df < 10.\n\n")
} else {
  cat(sprintf("NOTE: %d observations have df < 10:\n", n_small_df))
  small_df <- DT[DT$flag_small_df == 1,
                  c("newid", "obs", "z", "pcc", "sez", "df", "tstat")]
  small_df[, c("z", "pcc", "sez", "tstat")] <-
    round(small_df[, c("z", "pcc", "sez", "tstat")], 4)
  print(small_df, row.names = FALSE)
  cat("\n")
}

# ============================================================
# SUMMARY
# ============================================================

cat("==========================================================\n")
cat("SUMMARY OF FLAGS\n")
cat("==========================================================\n")

flag_vars <- c("flag_true_duplicate", "flag_coding_conflict",
               "flag_se_minor", "flag_se_major",
               "flag_outlier_z1", "flag_outlier_z2", "flag_small_df")

for (fv in flag_vars) {
  n     <- sum(DT[[fv]])
  pct   <- 100 * n / nrow(DT)
  label <- switch(fv,
    flag_true_duplicate  = "True duplicate (drop one copy)",
    flag_coding_conflict = "Coding conflict (same estimate, diff. coding)",
    flag_se_minor        = "SE deficit minor (0.001-0.010)",
    flag_se_major        = "SE deficit major (> 0.010)",
    flag_outlier_z1      = "Outlier: |z| > 1.0",
    flag_outlier_z2      = "Outlier: |z| > mean +/- 3 SD",
    flag_small_df        = "Very small df (< 10)"
  )
  cat(sprintf("  %-45s %4d obs  (%5.1f%%)\n", label, n, pct))
}

# Combined flag: any issue at all
DT$flag_any <- as.integer(
  DT$flag_true_duplicate  == 1 |
  DT$flag_coding_conflict == 1 |
  DT$flag_se_major        == 1 |
  DT$flag_outlier_z1      == 1 |
  DT$flag_small_df        == 1
)
# (Minor SE deficit is excluded from flag_any as it is not
# typically actionable -- see Check 6 discussion above.)

cat(sprintf("\n  %-45s %4d obs  (%5.1f%%)\n",
            "Any flag (excl. minor SE deficit)",
            sum(DT$flag_any), 100 * mean(DT$flag_any)))
cat("\nNote: Flags identify observations for review, not automatic deletion.\n")
cat("      The researcher should inspect flagged observations and decide\n")
cat("      whether to exclude, retain, or correct them.\n\n")

# ============================================================
# Export: diagnostic report as Excel and flagged data as .dta
#
# Five Excel sheets:
#   1. Summary -- counts and percentages for each flag
#   2. SE deficit -- all observations where se_deficit > 0.001
#   3. Outliers -- flagged effect sizes and small df
#   4. True duplicates -- duplicate groups for review
#   5. Coding conflicts -- conflicting rows and the variables that differ
#
# The full dataset with flag columns is also saved as a new
# .dta file for researcher inspection only. Analysis scripts
# do NOT load it; they load SCData_processed.rds (saved below).
# ============================================================

# --- Sheet 1: Summary ---

summary_sheet <- tibble(
  Flag = c(
    "True duplicate (identical on all vars except obs)",
    "Coding conflict (same z/sez/tstat, different coding)",
    "SE deficit minor (0.001 to 0.010)",
    "SE deficit major (> 0.010)",
    "Outlier: |z| > 1.0",
    "Outlier: |z| > mean +/- 3 SD",
    "Very small df (< 10)",
    "Any flag (excl. minor SE deficit)"
  ),
  N_obs = c(
    sum(DT$flag_true_duplicate),
    sum(DT$flag_coding_conflict),
    sum(DT$flag_se_minor),
    sum(DT$flag_se_major),
    sum(DT$flag_outlier_z1),
    sum(DT$flag_outlier_z2),
    sum(DT$flag_small_df),
    sum(DT$flag_any)
  ),
  Pct_of_total = round(c(
    mean(DT$flag_true_duplicate),
    mean(DT$flag_coding_conflict),
    mean(DT$flag_se_minor),
    mean(DT$flag_se_major),
    mean(DT$flag_outlier_z1),
    mean(DT$flag_outlier_z2),
    mean(DT$flag_small_df),
    mean(DT$flag_any)
  ) * 100, 1),
  N_studies = c(
    length(unique(DT$newid[DT$flag_true_duplicate  == 1])),
    length(unique(DT$newid[DT$flag_coding_conflict == 1])),
    length(unique(DT$newid[DT$flag_se_minor        == 1])),
    length(unique(DT$newid[DT$flag_se_major        == 1])),
    length(unique(DT$newid[DT$flag_outlier_z1      == 1])),
    length(unique(DT$newid[DT$flag_outlier_z2      == 1])),
    length(unique(DT$newid[DT$flag_small_df        == 1])),
    length(unique(DT$newid[DT$flag_any             == 1]))
  )
)

# --- Sheet 2: SE deficit detail ---

se_detail <- DT[DT$se_deficit > 0.001,
                 c("newid", "obs", "z", "pcc", "sez", "df",
                   "sez_theoretical", "se_deficit",
                   "flag_se_minor", "flag_se_major")] %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(newid, obs)

# --- Sheet 3: Outliers and small df ---

out_detail <- DT[DT$flag_outlier_z1 == 1 | DT$flag_outlier_z2 == 1 | DT$flag_small_df == 1,
                  c("newid", "obs", "z", "pcc", "sez", "df", "tstat",
                    "flag_outlier_z1", "flag_outlier_z2", "flag_small_df")] %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(newid, obs)

# --- Sheet 4: True duplicates ---
# Add a duplicate_group column so it is immediately clear which
# rows are copies of each other. Rows sharing the same group
# number within a study are identical and one should be dropped.

true_dup_detail <- DT %>%
  filter(flag_true_duplicate == 1) %>%
  group_by(across(all_of(all_vars_except_obs))) %>%
  mutate(duplicate_group = cur_group_id()) %>%
  ungroup() %>%
  select(duplicate_group, newid, obs, z, pcc, sez, df, tstat) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(duplicate_group, obs)

# --- Sheet 5: Coding conflicts ---
# For each conflicting group, show the obs identifiers and
# the variables that differ within the group.

conflict_detail <- DT[DT$flag_coding_conflict == 1, ] %>%
  select(newid, obs, z, sez, tstat,
         SC1_Cognitive, SC1_Structural, SC1_Other,
         SC2_GenTrust, SC2_PartTrust, SC2_AssocGroup,
         SC2_SocialPart, SC2_Bonding, SC2_Bridging,
         SC2_EthnicFrag, SC2_OtherSC, sc_category2,
         Endog_IV, Endog_FE, NoEndogeneity) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  arrange(newid, obs)

write_xlsx(
  list(
    "Summary"          = summary_sheet,
    "SE deficit"       = se_detail,
    "Outliers"         = out_detail,
    "True duplicates"  = true_dup_detail,
    "Coding conflicts" = conflict_detail
  ),
  "DataCleaning_Report.xlsx"
)

# --- Save full dataset with flags as a new .dta file ---
# This file is saved for researcher inspection only (e.g.,
# to review flagged observations in Stata). It is NOT loaded
# by any subsequent analysis scripts -- those scripts load
# SCData_processed.rds instead.

write_dta(DT, "SCData_flagged.dta")

# --- Save processed dataset ---
# For this protocol, the processed dataset is identical to the raw data:
# no flagged observations are dropped, in order to maintain comparability
# with Xue, Reed, and van Aert (2025). In your own work, inspect
# DataCleaning_Report.xlsx, decide which flagged observations to drop or
# revise, make those changes to DT, and then save the result here.
# Keep a log file documenting every change made and the reason for it.
saveRDS(DT, "SCData_processed.rds")

cat("==========================================================\n")
cat("Output saved:\n")
cat("  DataCleaning_Report.xlsx  -- diagnostic report\n")
cat("  SCData_flagged.dta        -- original data + flag columns\n")
cat("  SCData_processed.rds      -- dataset for all subsequent analysis\n")
cat("                               (identical to raw data in this protocol;\n")
cat("                                modify before saving for your own work)\n")
cat("==========================================================\n")

# ── Run time ──────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["el