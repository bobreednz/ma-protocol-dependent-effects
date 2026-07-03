library(haven)
library(RoBMA)
library(readxl)

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# ---- LOAD FITTED RoBMA MODEL FROM PROTOCOL_FIGURE3.R ----
# Avoids refitting the model (~5 hours). Must run Protocol_Figure3.R first.
cat("Loading intercept-only RoBMA model...\n")
fit <- readRDS("RoBMA_intercept_only.rds")

# ---- EXTRACT RoBMA BIAS-CORRECTED MEAN ----
# summary() returns model-averaged posterior estimates.
# The "mu" row gives the bias-corrected mean effect (Fisher's z scale).
s <- summary(fit)

# The estimates table has rows named by parameter; extract mu's posterior mean
# and 95% credible interval.
est_table  <- s$estimates
mu_row     <- est_table["mu", ]
robma_mean <- mu_row[["Mean"]]
robma_lo   <- mu_row[["0.025"]]   # lower bound of 95% CI
robma_hi   <- mu_row[["0.975"]]   # upper bound of 95% CI

cat(sprintf("RoBMA bias-corrected mean (Fisher's z): %.3f [%.3f, %.3f]\n",
            robma_mean, robma_lo, robma_hi))

# ---- READ CHE MEAN FROM TABLE 4 ----
# Table 4 has five estimators; we use CHE (the last column).
# Row 1 is the header, Row 2 is the Constant (intercept).
t4 <- read_excel("Table4_Protocol.xlsx", col_names = TRUE)

# The CHE column is the 6th column; the intercept is in the second data row.
# We strip significance stars and convert to numeric.
che_mean_raw <- t4[[6]][1]                          # "0.175***"
che_mean     <- as.numeric(gsub("[^0-9.-]", "", che_mean_raw))
che_se_raw   <- t4[[6]][2]                          # "(0.018)"
che_se       <- as.numeric(gsub("[^0-9.-]", "", che_se_raw))

cat(sprintf("CHE uncorrected mean (Fisher's z): %.3f (SE = %.3f)\n",
            che_mean, che_se))

# ---- READ FAT-PET CORRECTED ESTIMATE FROM TABLE 5 ----
# Table 5 FAT-PET results; "Effect beyond bias" row, CHE column.
# IMPORTANT: must read the "Panel B - Full controls" sheet explicitly. Table 5
# is saved with Panel A ("No controls") as the first sheet, so calling
# read_excel() without a sheet argument silently reads Panel A instead of the
# Panel B estimate that Section 6's decision rule designates as preferred
# (the CHE PET intercept in Panel B, not Panel A -- see the "PET vs PEESE"
# discussion in the protocol). This was previously a bug: the comparison
# table below was built from Panel A's 0.028 instead of Panel B's 0.042.
t5 <- read_excel("Table5_Protocol.xlsx", sheet = "Panel B - Full controls", col_names = TRUE)

fatpet_mean_raw <- t5[[6]][1]                       # "0.042"
fatpet_mean     <- as.numeric(gsub("[^0-9.-]", "", fatpet_mean_raw))
fatpet_se_raw   <- t5[[6]][2]                       # "(0.032)"
fatpet_se       <- as.numeric(gsub("[^0-9.-]", "", fatpet_se_raw))

cat(sprintf("FAT-PET corrected mean (Fisher's z): %.3f (SE = %.3f)\n",
            fatpet_mean, fatpet_se))

# ---- BUILD COMPARISON TABLE ----
# All estimates are on the Fisher's z scale.
# RoBMA reports a posterior mean with a 95% credible interval;
# CHE and FAT-PET report a point estimate with a standard error.

comparison <- data.frame(
  Estimator   = c("CHE (uncorrected)", "FAT-PET (CHE)", "RoBMA"),
  Mean_Fishers_z = c(che_mean, fatpet_mean, robma_mean),
  Uncertainty = c(
    sprintf("SE = %.3f", che_se),
    sprintf("SE = %.3f", fatpet_se),
    sprintf("95%% CI [%.3f, %.3f]", robma_lo, robma_hi)
  ),
  stringsAsFactors = FALSE
)

cat("\n=== Comparison Table (Fisher's z scale) ===\n")
print(comparison, row.names = FALSE)

# ---- DECISION RULE ----
# If the FAT-PET and RoBMA corrected means tell the same substantive story,
# stop here. If they diverge, proceed to Tables 9-11 (Protocol_Table9.R onward).
cat("\n--- Decision Rule ---\n")
cat(sprintf("CHE uncorrected mean:    %.3f\n", che_mean))
cat(sprintf("FAT-PET corrected mean:  %.3f\n", fatpet_mean))
cat(sprintf("RoBMA corrected mean:    %.3f\n", robma_mean))

divergence <- abs(fatpet_mean - robma_mean)
cat(sprintf("Divergence (FAT-PET vs RoBMA): %.3f\n", divergence))
cat("If FAT-PET and RoBMA corrected estimates diverge substantively,\n")
cat("proceed to Protocol_Table9.R (BIC-based moderator selection).\n")

# ---- EXTRACT BAYES FACTORS FROM RoBMA ----
# summary()$inclusion_components reports evidence for each model component:
# effect (mu != 0), heterogeneity (tau > 0), and publication bias.
# The inclusion Bayes factor (BF) quantifies evidence in favor of each component.
components <- s$inclusion_components

cat("\n=== RoBMA Component Bayes Factors ===\n")
print(components)

# Extract by column position to avoid name-matching issues.
# Column order from summary output: Prior prob., Post. prob., Inclusion BF, error%
# Row order: Effect, Heterogeneity, Publication Bias
pip_col <- 2   # Post. prob.
bf_col  <- 3   # Inclusion BF

bf_effect  <- components[1, bf_col]
bf_hetero  <- components[2, bf_col]
bf_pub     <- components[3, bf_col]

pip_effect <- components[1, pip_col]
pip_hetero <- components[2, pip_col]
pip_pub    <- components[3, pip_col]

bf_table <- data.frame(
  Component           = c("Effect (mu != 0)", "Heterogeneity (tau > 0)", "Publication bias"),
  Inclusion_BF        = c(bf_effect, bf_hetero, bf_pub),
  Posterior_Inclusion_Prob = c(pip_effect, pip_hetero, pip_pub),
  stringsAsFactors    = FALSE
)

cat("\n=== Bayes Factor Table ===\n")
print(bf_table, row.names = FALSE)

# ---- SAVE BOTH PANELS TO XLSX ----
library(openxlsx)

wb <- createWorkbook()

# Panel 1: mean effect comparison across estimators
addWorksheet(wb, "Comparison")
writeData(wb, "Comparison", comparison)

# Panel 2: RoBMA Bayes factors for effect, heterogeneity, publication bias
addWorksheet(wb, "BayesFactors")
writeData(wb, "BayesFactors", bf_table)

saveWorkbook(wb, "Table8_RoBMA_Comparison.xlsx", overwrite = TRUE)
cat("Table saved to Table8_RoBMA_Comparison.xlsx (two sheets)\