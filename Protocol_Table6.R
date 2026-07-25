# =============================================================================
# Protocol_Table6.R
# =============================================================================
# Purpose: Full meta-regression with all 19 regressors (sez + 18 study
# characteristics), combining Bayesian Model Averaging (BMA) and a full
# correlated-and-hierarchical effects (CHE) model. Replicates Table 10 in
# Xue et al. (2025), replacing the OLS full-model columns with CHE estimates.
#
# Table structure:
#   Panel A -- Estimates (7 columns):
#     Moderator | BMA PIP | BMA Post Mean | BMA Post SD |
#     CHE Coeff | CHE SE | CHE p-value
#     BMA columns are descriptive only (OLS-based, ignores clustering and
#     inverse-variance weighting). CHE columns are the basis for inference.
#
#   Panel B -- Interpretation (2 columns):
#     Moderator | Interpretation
#     Anchored on CHE significance (p < 0.10). BMA posterior dispersion
#     (SD > |mean|) is noted where relevant as a flag for specification
#     sensitivity. sez is excluded (methodological control, not a substantive
#     moderator). Sorted: Evidence of Moderator Effect first, No Evidence of
#     Moderator Effect below; alphabetical within each group.
#
# BMA approach:
#   BMA is estimated via ordinary least squares using the BMS package. This
#   is intentional: BMA with multilevel or meta-analytic variance structures
#   is not supported by available packages, so OLS-based BMA is the standard
#   approach in the meta-regression literature. The variable sez is forced
#   into every model via fixed.reg, since the publication bias proxy should
#   always be controlled for. All 18 study characteristics compete freely.
#
#   BMS settings:
#     g = "hyper=UIP"   -- unit information prior on regression coefficients
#     mprior = "random" -- Beta-Binomial(1,1) prior on model size
#     burn  = 50000     -- burn-in iterations (discarded)
#     iter  = 100000    -- post-burn-in iterations retained for inference
#     set.seed(12345)   -- for reproducibility
#
#   From the BMA output we extract, for each variable:
#     PIP       -- posterior inclusion probability
#     Post Mean -- posterior mean (averaged over all models, zero weight when
#                  variable excluded)
#     Post SD   -- posterior SD (descriptive measure of posterior dispersion
#                  across the BMA exercise; not a valid standard error given
#                  the misspecified likelihood)
#
# CHE approach:
#   The CHE model includes all 19 regressors simultaneously. The block-
#   diagonal covariance matrix V_mat is built with rho = 0.5. REML is used
#   for variance components. CR2 clustered standard errors are clustered at
#   the study level (newid). CHE coefficients, SEs, and p-values are the
#   primary basis for inference and for best-practice predictions in Table 7.
#
# Output:
#   Table6_Protocol.xlsx -- two sheets: Panel A and Panel B
# =============================================================================

# -- Packages -----------------------------------------------------------------

packages <- c("here", "metafor", "clubSandwich", "tidyverse", "openxlsx", "BMS")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)           # resolves file paths relative to the project root (.Rproj)
library(metafor)       # rma.mv() for CHE meta-regression
library(clubSandwich)  # coef_test() for CR2 standard errors
library(tidyverse)     # data manipulation
library(openxlsx)      # Excel export with formatting
library(BMS)           # bms() for Bayesian Model Averaging

# Record start time.
start_time <- proc.time()

# -- Data ---------------------------------------------------------------------

dat <- readRDS(here("SCData_processed.rds"))

cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# =============================================================================
# Define the 19 regressors
# =============================================================================

regressors <- c(
  "sez",            # publication bias proxy (forced into every BMA model)
  "SC1_Cognitive",  # social capital type: cognitive
  "SC1_Structural", # social capital type: structural (reference: SC1_Other)
  "DV_GrowthRate",  # outcome: GDP growth rate (reference: DV_GDPLevel)
  "PubYear",        # publication year (continuous)
  "Published",      # 1 = published in peer-reviewed journal
  "LaggedDV",       # 1 = lagged dependent variable included as control
  "LaggedSC",       # 1 = lagged social capital used
  "NumberSCVars",   # number of social capital variables in the regression
  "Endog_IV",       # 1 = instrumental variable used for endogeneity
  "Endog_FE",       # 1 = fixed effects used for endogeneity
  "CityLevel",      # 1 = city-level analysis
  "RegionLevel",    # 1 = region-level analysis
  "CountryLevel",   # 1 = country-level analysis (reference: OtherLevel)
  "PanelData",      # 1 = panel data (reference: cross-section)
  "Reg_OECDEurope", # 1 = OECD / European sample
  "Reg_US",         # 1 = US sample
  "Reg_Africa",     # 1 = African sample
  "Reg_Asia"        # 1 = Asian sample (reference: Reg_Other)
)

cat(sprintf("\n%d regressors: sez + 18 study characteristics.\n",
            length(regressors)))

# =============================================================================
# Part 1: Bayesian Model Averaging (BMS)
# =============================================================================

bma_dat <- dat %>%
  select(z, all_of(regressors)) %>%
  drop_na() %>%
  mutate(across(everything(), as.numeric))  # strip haven labels; bms() needs plain numeric

cat(sprintf("\nBMA dataset: %d complete observations.\n", nrow(bma_dat)))

set.seed(12345)

cat("Running BMA (burn = 50,000; iter = 100,000). This may take a few minutes.\n")

bma_fit <- bms(
  X.data    = bma_dat,
  burn      = 50000,
  iter      = 100000,
  g         = "hyper=UIP",
  mprior    = "random",
  fixed.reg = "sez",
  nmodel    = 2000,
  mcmc      = "bd"
)

bma_coef <- coef(bma_fit, exact = TRUE)

# Extract by column name -- robust to column ordering differences across BMS versions.
bma_post_mean <- bma_coef[, "Post Mean"]
bma_post_sd   <- bma_coef[, "Post SD"]
bma_pip       <- bma_coef[, "PIP"]

bma_table <- data.frame(
  variable  = rownames(bma_coef),
  post_mean = as.numeric(bma_post_mean),
  post_sd   = as.numeric(bma_post_sd),
  pip       = as.numeric(bma_pip)
)

cat("BMA complete.\n")
cat(sprintf("  Variables with PIP > 0.50: %d\n",
            sum(bma_table$pip > 0.50, na.rm = TRUE)))

# =============================================================================
# Part 2: CHE full meta-regression
# =============================================================================

# -- 2a. Block-diagonal covariance matrix (rho = 0.5) ------------------------

rho <- 0.5

V_mat <- vcalc(
  vi      = dat$sez^2,
  cluster = dat$newid,
  rho     = rho
)

# -- 2b. Build the moderator formula -----------------------------------------

mods_formula <- as.formula(
  paste("~", paste(regressors, collapse = " + "))
)

cat("\nFitting CHE full meta-regression ...\n")
cat(sprintf("  Formula: %s\n", deparse(mods_formula)))

# -- 2c. Fit CHE with REML ---------------------------------------------------

che_full <- tryCatch(
  rma.mv(
    yi     = z,
    V      = V_mat,
    mods   = mods_formula,
    random = ~ 1 | newid/obs,
    data   = dat,
    method = "REML"
  ),
  error = function(e) {
    cat(sprintf("ERROR fitting CHE: %s\n", conditionMessage(e)))
    NULL
  }
)

if (is.null(che_full)) stop("CHE model failed to converge. Check the data.")

cat("CHE model fitted successfully.\n")
cat(sprintf("  tau  = %.4f  (between-study SD)\n", sqrt(che_full$sigma2[1])))
cat(sprintf("  omega= %.4f  (within-study SD)\n",  sqrt(che_full$sigma2[2])))

# -- 2d. CR2 standard errors -------------------------------------------------

ct <- coef_test(
  obj     = che_full,
  vcov    = "CR2",
  cluster = dat$newid
)

ct_df <- as.data.frame(ct)
ct_df$variable_raw <- rownames(ct_df)

cat("\nCHE coefficients (CR2 SEs):\n")
print(ct_df[, c("variable_raw", "beta", "SE", "tstat", "p_Satt")])

# =============================================================================
# Part 3: Assemble Panel A and Panel B
# =============================================================================

# -- Nicer variable labels ----------------------------------------------------

var_labels <- c(
  sez            = "Standard error (sez)",
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

# -- Align BMA and CHE rows to regressors ordering ---------------------------

bma_ordered <- bma_table[match(regressors, bma_table$variable), ]
rownames(bma_ordered) <- NULL

# Remove intercept row from coef_test() output and strip any "mods" prefix.
ct_mods <- ct_df[ct_df$variable_raw != "intrcpt", ]
ct_mods$variable <- sub("^mods", "", ct_mods$variable_raw)
ct_ordered <- ct_mods[match(regressors, ct_mods$variable), ]
rownames(ct_ordered) <- NULL

# -- Helper functions ---------------------------------------------------------

# Round to 3 decimal places; NA becomes blank.
fmt3 <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "f"))

# Format p-value: show "< 0.001" for very small values.
fmt_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f")))
}

# Significance stars based on two-tailed Satterthwaite p-value.
stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.01  ~ "***",
    p < 0.05  ~ "**",
    p < 0.10  ~ "*",
    TRUE      ~ ""
  )
}

# -- Panel A ------------------------------------------------------------------
# All 19 moderators, original ordering.
# CHE coefficient has significance stars appended.

panel_a <- data.frame(
  Moderator      = var_labels[regressors],
  `BMA PIP`      = fmt3(bma_ordered$pip),
  `BMA Post Mean`= fmt3(bma_ordered$post_mean),
  `BMA Post SD`  = fmt3(bma_ordered$post_sd),
  `CHE Coeff`    = paste0(fmt3(ct_ordered$beta), stars(ct_ordered$p_Satt)),
  `CHE SE`       = paste0("(", fmt3(ct_ordered$SE), ")"),
  `CHE p-value`  = fmt_p(ct_ordered$p_Satt),
  check.names    = FALSE,
  stringsAsFactors = FALSE
)

# -- Panel B ------------------------------------------------------------------
# 18 moderators (sez excluded). Interpretation anchored on CHE significance
# (p < 0.10). BMA posterior dispersion is flagged where PIP >= 0.10 and
# SD > |Post Mean|, indicating the effect size estimate varies substantially
# across model specifications.

moderators_only <- regressors[regressors != "sez"]

panel_b_raw <- data.frame(
  variable   = moderators_only,
  label      = var_labels[moderators_only],
  che_p      = ct_ordered$p_Satt[match(moderators_only, regressors)],
  bma_mean   = bma_ordered$post_mean[match(moderators_only, regressors)],
  bma_sd     = bma_ordered$post_sd[match(moderators_only, regressors)],
  bma_pip    = bma_ordered$pip[match(moderators_only, regressors)],
  stringsAsFactors = FALSE
)

# Evidence of Moderator Effect = CHE p < 0.10.
panel_b_raw$supported <- !is.na(panel_b_raw$che_p) & panel_b_raw$che_p < 0.10

# High-dispersion flag: PIP >= 0.10 and BMA posterior SD > |posterior mean|.
# Applied only when posterior mean is non-negligible (avoids spurious flags
# when BMA has essentially zeroed out a variable).
panel_b_raw$high_disp <- with(panel_b_raw,
  !is.na(bma_pip) & bma_pip >= 0.10 &
  !is.na(bma_sd)  & !is.na(bma_mean) & abs(bma_mean) > 1e-6 &
  bma_sd / abs(bma_mean) > 1
)

# Build interpretation text.
panel_b_raw$interpretation <- with(panel_b_raw, {
  base <- ifelse(supported, "Evidence of Moderator Effect", "No Evidence of Moderator Effect")
  disp_note <- ifelse(high_disp,
    "; BMA posterior SD > abs(mean) (effect sensitive to model specification)",
    "")
  paste0(base, disp_note)
})

# Sort: Evidence of Moderator Effect first, then No Evidence of Moderator Effect; alphabetical by label within groups.
panel_b_sorted <- panel_b_raw[order(!panel_b_raw$supported, panel_b_raw$label), ]

panel_b <- data.frame(
  Moderator      = panel_b_sorted$label,
  Interpretation = panel_b_sorted$interpretation,
  stringsAsFactors = FALSE
)

# -- Console summaries --------------------------------------------------------

cat("\n-- BMA summary ----------------------------------------------------------\n")
cat("Variables with PIP >= 0.50:\n")
high_pip <- bma_table[bma_table$pip >= 0.50, ]
if (nrow(high_pip) == 0) {
  cat("  None.\n")
} else {
  for (j in seq_len(nrow(high_pip))) {
    cat(sprintf("  %-20s  PIP = %.3f  PostMean = %.4f  PostSD = %.4f\n",
                high_pip$variable[j], high_pip$pip[j],
                high_pip$post_mean[j], high_pip$post_sd[j]))
  }
}

cat("\n-- CHE summary ----------------------------------------------------------\n")
cat("Coefficients significant at 10% level:\n")
sig_che <- ct_mods[!is.na(ct_mods$p_Satt) & ct_mods$p_Satt < 0.10, ]
if (nrow(sig_che) == 0) {
  cat("  None.\n")
} else {
  for (j in seq_len(nrow(sig_che))) {
    cat(sprintf("  %-20s  beta = %.4f  SE = %.4f  p = %.4f\n",
                sig_che$variable[j], sig_che$beta[j],
                sig_che$SE[j], sig_che$p_Satt[j]))
  }
}

cat("\n-- Panel B: Interpretation ----------------------------------------------\n")
print(panel_b)

# =============================================================================
# Export to Excel (openxlsx) -- two sheets: Panel A and Panel B
# =============================================================================

wb <- createWorkbook()
addWorksheet(wb, "Panel A")
addWorksheet(wb, "Panel B")

# -- Shared styles ------------------------------------------------------------

style_title       <- createStyle(fontSize = 11, textDecoration = "bold")
style_header      <- createStyle(fontSize = 10, textDecoration = "bold",
                                  border = "Bottom", borderStyle = "medium",
                                  wrapText = TRUE, halign = "center")
style_header_left <- createStyle(fontSize = 10, textDecoration = "bold",
                                  border = "Bottom", borderStyle = "medium",
                                  wrapText = TRUE, halign = "left")
style_group       <- createStyle(fontSize = 10, textDecoration = "italic",
                                  halign = "center")
style_note        <- createStyle(fontSize = 9, wrapText = TRUE)
style_center      <- createStyle(halign = "center")

# =============================================================================
# Sheet 1: Panel A -- Estimates
# =============================================================================

r <- 1  # current row tracker for Panel A

# Title
writeData(wb, "Panel A", "Panel A: Meta-Regression Estimates",
          startRow = r, startCol = 1)
mergeCells(wb, "Panel A", rows = r, cols = 1:7)
addStyle(wb, "Panel A", style_title, rows = r, cols = 1)
r <- r + 1

# Group header: write each non-empty cell individually to avoid NA from
# empty-string data.frame coercion.
writeData(wb, "Panel A", "BMA (descriptive only)", startRow = r, startCol = 2)
writeData(wb, "Panel A", "CHE (primary)",          startRow = r, startCol = 5)
mergeCells(wb, "Panel A", rows = r, cols = 2:4)
mergeCells(wb, "Panel A", rows = r, cols = 5:7)
addStyle(wb, "Panel A", style_group, rows = r, cols = 2:7, gridExpand = TRUE)
r <- r + 1

# Column headers
writeData(wb, "Panel A", "Moderator",        startRow = r, startCol = 1)
writeData(wb, "Panel A", "BMA: PIP",        startRow = r, startCol = 2)
writeData(wb, "Panel A", "BMA: Mean",       startRow = r, startCol = 3)
writeData(wb, "Panel A", "BMA: SD",         startRow = r, startCol = 4)
writeData(wb, "Panel A", "CHE: Coeff",      startRow = r, startCol = 5)
writeData(wb, "Panel A", "CHE: SE",         startRow = r, startCol = 6)
writeData(wb, "Panel A", "CHE: p-value",    startRow = r, startCol = 7)
addStyle(wb, "Panel A", style_header_left, rows = r, cols = 1)
addStyle(wb, "Panel A", style_header,      rows = r, cols = 2:7, gridExpand = TRUE)
r <- r + 1

# Data rows
writeData(wb, "Panel A", panel_a, startRow = r, startCol = 1, colNames = FALSE)
addStyle(wb, "Panel A", style_center,
         rows = r:(r + nrow(panel_a) - 1), cols = 2:7, gridExpand = TRUE)
r <- r + nrow(panel_a) + 1  # +1 for blank separator

# Footer (CHE model statistics)
n_obs     <- nrow(dat)
n_studies <- length(unique(dat$newid))
tau_val   <- sqrt(che_full$sigma2[1])
omega_val <- sqrt(che_full$sigma2[2])

writeData(wb, "Panel A", "Observations",   startRow = r,     startCol = 1)
writeData(wb, "Panel A", n_obs,            startRow = r,     startCol = 5)
writeData(wb, "Panel A", "Studies",        startRow = r + 1, startCol = 1)
writeData(wb, "Panel A", n_studies,        startRow = r + 1, startCol = 5)
writeData(wb, "Panel A", "rho (assumed)",  startRow = r + 2, startCol = 1)
writeData(wb, "Panel A", fmt3(rho),        startRow = r + 2, startCol = 5)
writeData(wb, "Panel A", "tau",            startRow = r + 3, startCol = 1)
writeData(wb, "Panel A", fmt3(tau_val),    startRow = r + 3, startCol = 5)
writeData(wb, "Panel A", "omega",          startRow = r + 4, startCol = 1)
writeData(wb, "Panel A", fmt3(omega_val),  startRow = r + 4, startCol = 5)
r <- r + 6  # 5 footer rows + 1 blank

# Notes
notes_a <- c(
  "*** p < 0.01, ** p < 0.05, * p < 0.10. CHE standard errors are CR2, clustered at the study level (Satterthwaite df).",
  "BMA: BMS package, g = hyper-UIP, mprior = random, burn = 50,000, iter = 100,000, seed = 12345. sez forced into every model.",
  "BMA columns are descriptive only. The BMS implementation uses OLS and does not account for inverse-variance weighting or within-study dependence.",
  "BMA informs interpretation; CHE determines inference."
)

for (note in notes_a) {
  writeData(wb, "Panel A", note, startRow = r, startCol = 1)
  mergeCells(wb, "Panel A", rows = r, cols = 1:7)
  addStyle(wb, "Panel A", style_note, rows = r, cols = 1)
  r <- r + 1
}

# Column widths
setColWidths(wb, "Panel A", cols = 1:7,
             widths = c(28, 8, 10, 10, 12, 10, 10))

# =============================================================================
# Sheet 2: Panel B -- Interpretation
# =============================================================================

r2 <- 1  # current row tracker for Panel B

# Title
writeData(wb, "Panel B",
          "Panel B: Interpretation (sez excluded; CHE p < 0.10 = Evidence of Moderator Effect)",
          startRow = r2, startCol = 1)
mergeCells(wb, "Panel B", rows = r2, cols = 1:2)
addStyle(wb, "Panel B", style_title, rows = r2, cols = 1)
r2 <- r2 + 1

# Column headers
writeData(wb, "Panel B", "Moderator",      startRow = r2, startCol = 1)
writeData(wb, "Panel B", "Interpretation", startRow = r2, startCol = 2)
addStyle(wb, "Panel B", style_header_left, rows = r2, cols = 1:2, gridExpand = TRUE)
r2 <- r2 + 1

# Data rows
writeData(wb, "Panel B", panel_b, startRow = r2, startCol = 1, colNames = FALSE)
r2 <- r2 + nrow(panel_b) + 1

# Note
writeData(wb, "Panel B",
          "Interpretation is anchored on CHE significance (p < 0.10). BMA posterior SD > |mean| (where PIP >= 0.10) is noted as an indicator of sensitivity to model specification.",
          startRow = r2, startCol = 1)
mergeCells(wb, "Panel B", rows = r2, cols = 1:2)
addStyle(wb, "Panel B", style_note, rows = r2, cols = 1)

# Column widths
setColWidths(wb, "Panel B", cols = 1:2, widths = c(28, 60))

# =============================================================================
# Save
# =============================================================================

saveWorkbook(wb, here("Table6_Protocol.xlsx"), overwrite = TRUE)

cat("\nTable6_Protocol.xlsx saved (two sheets: Panel A, Panel B).\n")

# -- Run time -----------------------------------------------------------------

elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds.\n", elapsed["elapsed"]))
