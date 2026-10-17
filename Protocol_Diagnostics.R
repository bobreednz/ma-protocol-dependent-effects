# =============================================================================
# Protocol_Diagnostics.R
# =============================================================================
# Purpose: two comparison tables. Sheet "SelectionWeights" is the Section 6.2
# table used to diagnose selection. Sheet "InterceptOnly" is not displayed in
# the protocol; Protocol_Sensitivity.R reads it for the Section 8 table. No
# models are fitted; the script only reads output saved by earlier scripts.
#
#   Sheet "SelectionWeights": the selection weights estimated by 3PSM, the
#       4PSM and RoBMA, lined up by one-sided p-value interval. Both measure how
#       likely an estimate in each interval is to be reported, relative to a
#       significant estimate of the expected sign.
#   Sheet "InterceptOnly": the uncorrected mean and every corrected mean
#       from Section 5 (before moderators are added), with the percentage reduction relative
#       to the uncorrected CHE mean.
#
# Inputs (run these scripts first):
#   Table_FiveEstimators.xlsx            Protocol_FiveEstimators.R    (uncorrected CHE mean)
#   Table_CHE_ISCW.xlsx             Protocol_CHE_ISCW.R  (CHE-PET-PEESE,
#                                                         CHE-ISCW-PET-PEESE)
#   Table_SelectionModels.xlsx      Protocol_SelectionModels.R
#   Table_RoBMA_CHE_Intercept.xlsx  Protocol_RoBMA_CHE_Intercept.R
#
# Output: Table_Diagnostics.xlsx
# Run time: a few seconds.
# =============================================================================

library(here)     # resolves file paths relative to the project root (.Rproj)
library(readxl)   # read_excel()
library(writexl)  # write_xlsx()

# ---- Helpers -----------------------------------------------------------------
# Table 4 stores formatted text such as "0.175***" and "[0.139, 0.212]".
# These helpers turn that text back into numbers.
num <- function(s) as.numeric(gsub("[^0-9.eE-]", "", s))
ci_parse <- function(s) as.numeric(strsplit(gsub("\\[|\\]|\\s", "", s), ",")[[1]])
cell <- function(tab, row, col = "CHE") tab[[col]][tab[[1]] == row][1]

f3  <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)
int <- function(lo, hi) if (is.na(lo)) "--" else sprintf("[%.3f, %.3f]", lo, hi)
red <- function(x, base) if (is.na(x)) "--" else sprintf("%.0f%%", 100 * (base - x) / base)

# ---- Read inputs ----------------------------------------------------------------
t4   <- read_excel(here("Table_FiveEstimators.xlsx"), col_types = "text")
iscw <- read_excel(here("Table_CHE_ISCW.xlsx"),          sheet = "Values")
selm <- read_excel(here("Table_SelectionModels.xlsx"),   sheet = "Values")
rob  <- read_excel(here("Table_RoBMA_CHE_Intercept.xlsx"), sheet = "Values")

getv <- function(tab, key, col = "Quantity") tab[tab[[col]] == key, ][1, ]

# ---- Uncorrected CHE mean (five-estimators table) ----------------------------
che     <- num(cell(t4, "Constant"))
che_ci  <- ci_parse(cell(t4, "95% CI"))

# PET-PEESE corrections, before moderators are added (the PET-PEESE rule has
# already been applied in Protocol_CHE_ISCW.R)
# Older versions of Table_CHE_ISCW.xlsx use the label "CHE PET-PEESE"
che_pp  <- getv(iscw, "CHE-PET-PEESE")
if (is.na(che_pp$Estimate)) che_pp <- getv(iscw, "CHE PET-PEESE")
iscw_pp <- getv(iscw, "CHE-ISCW-PET-PEESE")
bma     <- getv(rob, "BMA mu")
rob_ma  <- getv(rob, "RoBMA mu, model-averaged")
rob_co  <- getv(rob, "RoBMA mu, conditional")
psm3    <- selm[selm$Model == "3PSM" & selm$Parameter == "beta", ][1, ]
psm4    <- selm[selm$Model == "4PSM" & selm$Parameter == "beta", ][1, ]

# ---- Sheet 1: corrected means before moderators are added ------------------
row_io <- function(approach, estimator, est, lo, hi, type, base = che, reduction = TRUE) {
  data.frame(Approach = approach, Estimator = estimator,
             Mean = f3(est), `95% interval` = int(lo, hi), `Interval type` = type,
             `Reduction from uncorrected CHE mean` = if (reduction) red(est, base) else "--",
             check.names = FALSE, stringsAsFactors = FALSE)
}

intercept_only <- rbind(
  row_io("No correction", "CHE (frequentist)", che, che_ci[1], che_ci[2], "CR2 CI", reduction = FALSE),
  row_io("No correction", "BMA.mv (Bayesian CHE)", bma$Estimate, bma$Lower, bma$Upper, "CrI", reduction = FALSE),
  row_io("PET-PEESE", paste0("CHE-PET-PEESE (", che_pp$Correction, " selected)"),
         che_pp$Estimate, che_pp$Lower, che_pp$Upper, "CR2 CI"),
  row_io("PET-PEESE", paste0("CHE-ISCW-PET-PEESE (robustness check, ", iscw_pp$Correction, " selected)"),
         iscw_pp$Estimate, iscw_pp$Lower, iscw_pp$Upper, "CR2 CI"),
  row_io("Selection model", "4PSM", psm4$Estimate, psm4$Lower, psm4$Upper, "Bootstrap CI"),
  row_io("Selection model", "3PSM", psm3$Estimate, psm3$Lower, psm3$Upper, "Bootstrap CI"),
  row_io("RoBMA", "Model-averaged", rob_ma$Estimate, rob_ma$Lower, rob_ma$Upper, "CrI"),
  row_io("RoBMA", "Conditional on an effect", rob_co$Estimate, rob_co$Lower, rob_co$Upper, "CrI")
)

# ---- Sheet 2: selection weights -----------------------------------------------
l1 <- selm[selm$Model == "4PSM" & selm$Parameter == "lambda1", ][1, ]
l3 <- selm[selm$Model == "3PSM" & selm$Parameter == "lambda1", ][1, ]   # 3PSM has one weight
l2 <- selm[selm$Model == "4PSM" & selm$Parameter == "lambda2", ][1, ]
w  <- function(q) { r <- getv(rob, q); paste(f3(r$Estimate), int(r$Lower, r$Upper)) }

# RoBMA splits the wrong-sign region (p > 0.5) into three intervals. When the
# three weights are the same (as here), they are reported once.
upper_q <- c("RoBMA omega[0.5,0.95]", "RoBMA omega[0.95,0.975]", "RoBMA omega[0.975,1]")
upper_e <- sapply(upper_q, function(q) getv(rob, q)$Estimate)
rob_wrong <- if (max(upper_e) - min(upper_e) < 0.001) {
  w(upper_q[1])
} else {
  paste0(w(upper_q[1]), " (0.5 to 0.95); ", w(upper_q[2]), " (0.95 to 0.975); ",
         w(upper_q[3]), " (0.975 to 1)")
}

# RoBMA has a threshold at 0.05 that the step-function models lack, so the
# expected-sign, not-significant range gets two rows. The 3PSM and 4PSM each
# estimate one weight for both rows ("Same as above").
selection_weights <- data.frame(
  `One-sided p-value` = c("0 to 0.025", "0.025 to 0.05", "0.05 to 0.5", "Above 0.5"),
  `What these estimates are` = c("Significant at 5% (two-sided), expected sign",
                                 "Expected sign, significant at 10% but not 5%",
                                 "Expected sign, not significant at 10%",
                                 "Wrong sign"),
  `3PSM weight [95% bootstrap CI]` = c("1 = reference",
                                       paste(f3(l3$Estimate), int(l3$Lower, l3$Upper)),
                                       "Same as above",
                                       "Same as above"),
  `4PSM weight [95% bootstrap CI]` = c("1 = reference",
                                       paste(f3(l1$Estimate), int(l1$Lower, l1$Upper)),
                                       "Same as above",
                                       paste(f3(l2$Estimate), int(l2$Lower, l2$Upper))),
  `RoBMA weight [95% CrI]` = c("1 = reference",
                               w("RoBMA omega[0.025,0.05]"),
                               w("RoBMA omega[0.05,0.5]"),
                               rob_wrong),
  check.names = FALSE, stringsAsFactors = FALSE
)

write_xlsx(list(SelectionWeights = selection_weights, InterceptOnly = intercept_only),
           path = here("Table_Diagnostics.xlsx"))
print(selection_weights, row.names = FALSE)
print(intercept_only, row.names = FALSE)
cat("\nSaved Table_Diagnostics.xlsx\n")
