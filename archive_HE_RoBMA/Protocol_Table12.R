# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Purpose
# ============================================================
#
# Table 12 pulls together, in one place, the key numbers that
# Section 10 uses to compare FAT-PET and RoBMA. It does not run
# any new models -- it simply reads the already-produced output
# files from Tables 4, 5, 7, 10, and 11 and assembles the seven
# numbers discussed in "The core comparison" into a single table.
# This keeps Section 10 from having to repeat those numbers in
# prose, and keeps the summary table automatically in sync with
# its source tables if any of them are ever re-run.
#
# Every row reports a 95% interval for the mean, rather than
# mixing standard errors, confidence intervals, and prediction
# intervals as an earlier version of this table did. A prediction
# interval answers a different question (the range a single
# future study's true effect might fall in, which is necessarily
# wider because it also incorporates heterogeneity) than a
# confidence or credible interval does (uncertainty about the
# mean itself), so mixing them in one column risked the reader
# misreading interval width as directly comparable across rows
# when it wasn't. Prediction intervals remain available in
# Tables 7 and 11 themselves.
#
# The FAT-PET rows report a 95% confidence interval (frequentist,
# built from the Satterthwaite df used throughout this protocol).
# The RoBMA rows report a 95% credible interval instead -- RoBMA
# is a Bayesian method, so its interval reflects the posterior
# distribution, not repeated-sampling coverage. The two are
# conceptually different even though both are often loosely
# called "confidence intervals" in applied writing, so the
# Interval Type column and the note below the table label them
# separately rather than calling everything "CI".
#
# ============================================================
# Install and load required packages
# ============================================================

packages <- c("here", "readxl", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)       # resolves file paths relative to the project root (.Rproj)
library(readxl)    # read the source Excel files
library(tidyverse) # data manipulation
library(writexl)   # export to Excel

# ============================================================
# 1. Uncorrected CHE mean (Table 4)
# ============================================================

tbl4 <- read_excel(here("Table4_Protocol.xlsx"), col_types = "text")

che_raw_est <- tbl4$CHE[tbl4$Variable == "Constant"]
che_raw_ci  <- tbl4$CHE[tbl4$Variable == "95% CI"]

# ============================================================
# 2. FAT-PET corrected mean, unconditional (Table 5, Panel B)
# ============================================================

tbl5b <- read_excel(here("Table5_Protocol.xlsx"), sheet = "Panel B - Full controls", col_types = "text")

fatpet_uncond_est <- tbl5b$CHE[tbl5b$Variable == "Effect beyond bias"]
fatpet_uncond_ci  <- tbl5b$CHE[tbl5b$Variable == "95% CI"]

# ============================================================
# 3. RoBMA corrected mean, unconditional (Table 10)
# ============================================================
#
# This is the average-study bias-corrected mean from the RoBMA
# meta-regression: pooled_effect(), the model-averaged (publication-
# bias-corrected) effect averaged across the moderators' actual
# empirical distribution in this dataset. As of 2026-07-27, Table 10
# fits both factor moderators with treatment contrasts and computes
# this quantity via pooled_effect() rather than marginal_means() under
# meandif contrasts, per Frantisek Bartos's direct advice that
# pooled_effect() -- not marginal_means() -- is the quantity
# corresponding to Table 5's CHE PET intercept, regardless of contrast
# coding. It is the RoBMA analogue of the FAT-PET Panel B intercept
# (both are with-controls corrected means at the sample's actual
# covariate composition), and it replaces the intercept-only Table 8
# figure used previously. Table 10 stores it (Fisher's z) with a 95%
# credible interval on its "AverageStudy" sheet.
# ============================================================

tbl10_avg <- read_excel(here("Table10_RoBMA_MetaRegression.xlsx"), sheet = "AverageStudy")

robma_uncond_est <- formatC(tbl10_avg$Mean[1], digits = 3, format = "f")
robma_uncond_ci  <- paste0("[", formatC(tbl10_avg$`2.5%`[1],  digits = 3, format = "f"),
                           ", ", formatC(tbl10_avg$`97.5%`[1], digits = 3, format = "f"), "]")

# ============================================================
# 4. Best-practice predictions: FAT-PET (Table 7) and RoBMA (Table 11)
#
# Both tables store one sheet per scenario (BP1, BP2), each with
# a "Fisher's z" row holding the Mean Prediction, 95% CI, and
# 95% PI. We use the 95% CI here (not the 95% PI) so that every
# row in this table reports the same kind of interval -- see the
# note at the top of this script.
# ============================================================

get_bp_row <- function(file, sheet) {
  tbl <- read_excel(here(file), sheet = sheet)
  tbl[tbl$`Effect Size` == "Fisher's z", ]
}

fatpet_bp1 <- get_bp_row("Table7_Protocol.xlsx", "BP1 - Non-OECD Europe")
fatpet_bp2 <- get_bp_row("Table7_Protocol.xlsx", "BP2 - OECD Europe")
robma_bp1  <- get_bp_row("Table11_RoBMA_BestPractice.xlsx", "BP1 - Non-OECD Europe")
robma_bp2  <- get_bp_row("Table11_RoBMA_BestPractice.xlsx", "BP2 - OECD Europe")

# ============================================================
# Assemble the summary table
# ============================================================

table_out <- tibble(
  Quantity = c(
    "Uncorrected mean effect",
    "Corrected mean effect (unconditional)",
    "Corrected mean effect (unconditional)",
    "Best-practice prediction: Non-OECD/Europe",
    "Best-practice prediction: Non-OECD/Europe",
    "Best-practice prediction: OECD/Europe",
    "Best-practice prediction: OECD/Europe"
  ),
  Method = c(
    "CHE (raw)", "FAT-PET", "RoBMA", "FAT-PET", "RoBMA", "FAT-PET", "RoBMA"
  ),
  `Fisher's z` = c(
    che_raw_est, fatpet_uncond_est, robma_uncond_est,
    fatpet_bp1$`Mean Prediction`, robma_bp1$`Mean Prediction`,
    fatpet_bp2$`Mean Prediction`, robma_bp2$`Mean Prediction`
  ),
  `Interval Type` = c(
    "95% CI", "95% CI", "95% CrI", "95% CI", "95% CrI", "95% CI", "95% CrI"
  ),
  `95% Interval` = c(
    che_raw_ci, fatpet_uncond_ci, robma_uncond_ci,
    fatpet_bp1$`95% CI`, robma_bp1$`95% CI`,
    fatpet_bp2$`95% CI`, robma_bp2$`95% CI`
  ),
  Source = c(
    "Table 4", "Table 5 (Panel B)", "Table 10",
    "Table 7", "Table 11", "Table 7", "Table 11"
  )
)

# ============================================================
# Export to Excel
# ============================================================

write_xlsx(table_out, here("Table12_ComparisonSummary.xlsx"))

cat("Done. Output saved to Table12_ComparisonSummary.xlsx\n")

# -- Run time -----------------------------------------------------------------
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds.\n", elapsed["elapsed"]))
