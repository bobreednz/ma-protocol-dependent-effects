# =============================================================================
# Protocol_Sensitivity.R
# =============================================================================
# Purpose: the summary table of Section 8. It collects, for each of the three
# corrections, the estimates reported in Sections 5 and 7:
#   - the corrected mean before moderators are added (Section 5);
#   - the corrected mean at average characteristics (Section 7);
#   - the two best-practice predictions, with prediction intervals (Section 7).
# The columns follow the fixed order of Sections 5 and 7; which correction is
# primary is stated in the text, so the table reads the same way whichever
# correction a researcher's diagnostics favor.
#
# The script fits no models. It reads the output files written by the
# scripts of Sections 5 to 7, so run those first:
#   Table_Diagnostics.xlsx          (Protocol_Diagnostics.R, sheet InterceptOnly)
#   Table_CHE_PETPEESE_MetaReg.xlsx (Protocol_CHE_PETPEESE_MetaReg.R; Section 7.2)
#   Table_CHE_PETPEESE_BestPractice.xlsx            (Protocol_CHE_PETPEESE_BestPractice.R; Section 7.2)
#   Table_Selection_MetaReg.xlsx    (Protocol_Selection_MetaReg.R; Section 7.3)
#   Table_RoBMA_CHE_MetaReg.xlsx    (Protocol_RoBMA_CHE_MetaReg.R; Section 7.4)
# The selection model (3PSM or 4PSM) is the one used in Section 7.3.
#
# Output: Table_Sensitivity.xlsx
#   "Summary"  the formatted table
#   "Values"   the same numbers, unformatted
# Run time: a few seconds.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
packages <- c("here", "readxl", "writexl")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)      # resolves file paths relative to the project root (.Rproj)
library(readxl)    # read_excel()
library(writexl)   # write_xlsx()

# ---- Helpers -------------------------------------------------------------------
f3  <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
iv  <- function(lo, hi) ifelse(is.na(lo), "--", sprintf("[%.3f, %.3f]", lo, hi))

# Turn a formatted cell such as "-0.022***" or "[-0.182, 0.139]" back into
# numbers (used for files that store formatted text only)
num <- function(txt) as.numeric(regmatches(txt, gregexpr("-?[0-9]+\\.[0-9]+", txt))[[1]])

# ---- Corrections before moderators (Section 5, via Table_Diagnostics.xlsx) ---
io <- read_excel(here("Table_Diagnostics.xlsx"), sheet = "InterceptOnly")
io_row <- function(pattern, approach = NULL) {
  keep <- grepl(pattern, io$Estimator)
  if (!is.null(approach)) keep <- keep & io$Approach == approach
  r <- io[keep, ]
  if (nrow(r) == 0) stop("Row not found in Table_Diagnostics.xlsx: ", pattern)
  r <- r[1, ]
  ci <- num(r$`95% interval`)
  c(est = as.numeric(r$Mean), lo = ci[1], hi = ci[2])
}

# ---- Selection model used in Section 7.3 ---------------------------------------
sel_vals <- read_excel(here("Table_Selection_MetaReg.xlsx"), sheet = "Values")
psm      <- sel_vals$Model[1]                       # "3PSM" or "4PSM"
sel_row  <- function(q) {
  r <- sel_vals[sel_vals$Quantity == q, ]
  c(est = r$Estimate, lo = r$Lower, hi = r$Upper, pi_lo = r$PI_lower, pi_hi = r$PI_upper)
}

# ---- CHE-PET-PEESE (Section 7.2) ------------------------------------------------
pet_vals <- read_excel(here("Table_CHE_PETPEESE_MetaReg.xlsx"), sheet = "Values")
pet_avg  <- pet_vals[pet_vals$Quantity == "Mean effect" & pet_vals$Selected, ]
# Table_CHE_PETPEESE_BestPractice.xlsx stores formatted text; read the Fisher's z row of each panel
pet_bp <- function(sheet) {
  r  <- read_excel(here("Table_CHE_PETPEESE_BestPractice.xlsx"), sheet = sheet)[1, ]
  ci <- num(r$`95% CI`); pi <- num(r$`95% PI`)
  c(est = num(r$`Mean Prediction`)[1], lo = ci[1], hi = ci[2], pi_lo = pi[1], pi_hi = pi[2])
}

# ---- RoBMA (Section 7.4) ---------------------------------------------------------
rob_vals <- read_excel(here("Table_RoBMA_CHE_MetaReg.xlsx"), sheet = "Values")
rob_row  <- function(q) {
  r <- rob_vals[rob_vals$Quantity == q, ]
  c(est = r$Estimate, lo = r$Lower, hi = r$Upper, pi_lo = r$PI_lower, pi_hi = r$PI_upper)
}

# ---- Assemble: one list per correction ---------------------------------------------
corr <- list(
  "CHE-PET-PEESE" = list(
    before = io_row("^CHE[- ]PET-PEESE"),
    avg    = c(est = pet_avg$Estimate, lo = pet_avg$Lower, hi = pet_avg$Upper),
    bp1    = pet_bp("BP1 - Other Regions"),
    bp2    = pet_bp("BP2 - OECD Europe")),
  "Selection model" = list(
    before = io_row(paste0("^", psm, "$")),
    avg    = sel_row("Corrected mean"),
    bp1    = sel_row("BP1 (Other Regions)"),
    bp2    = sel_row("BP2 (OECD/Europe)")),
  "RoBMA" = list(
    before = io_row("^Model-averaged$", approach = "RoBMA"),
    avg    = rob_row("RoBMA mean at average characteristics"),
    bp1    = rob_row("RoBMA BP1 (Other Regions)"),
    bp2    = rob_row("RoBMA BP2 (OECD/Europe)"))
)

# Column headings: the selection model is named (3PSM or 4PSM)
headings <- c("CHE-PET-PEESE", psm, "RoBMA")

# ---- Sheet 1: the formatted table --------------------------------------------------
column <- function(x) c(
  f3(x$before["est"]), iv(x$before["lo"], x$before["hi"]),
  f3(x$avg["est"]),    iv(x$avg["lo"], x$avg["hi"]),
  f3(x$bp1["est"]),    iv(x$bp1["lo"], x$bp1["hi"]), iv(x$bp1["pi_lo"], x$bp1["pi_hi"]),
  f3(x$bp2["est"]),    iv(x$bp2["lo"], x$bp2["hi"]), iv(x$bp2["pi_lo"], x$bp2["pi_hi"]))

summary_tbl <- data.frame(
  Estimate = c("Mean before moderators (Section 5)", "95% interval",
               "Mean at average characteristics (Section 7)", "95% interval",
               "Best practice: Other Regions", "95% interval", "95% prediction interval",
               "Best practice: OECD/Europe", "95% interval", "95% prediction interval"),
  stringsAsFactors = FALSE)
for (i in seq_along(corr)) summary_tbl[[headings[i]]] <- column(corr[[i]])

# ---- Sheet 2: unformatted values ------------------------------------------------------
values <- do.call(rbind, lapply(seq_along(corr), function(i) {
  x <- corr[[i]]
  do.call(rbind, lapply(c("before", "avg", "bp1", "bp2"), function(q) {
    v <- x[[q]]
    data.frame(Correction = headings[i], Quantity = q,
               Estimate = unname(v["est"]), Lower = unname(v["lo"]), Upper = unname(v["hi"]),
               PI_lower = if ("pi_lo" %in% names(v)) unname(v["pi_lo"]) else NA,
               PI_upper = if ("pi_hi" %in% names(v)) unname(v["pi_hi"]) else NA,
               stringsAsFactors = FALSE)
  }))
}))

write_xlsx(list(Summary = summary_tbl, Values = values), path = here("Table_Sensitivity.xlsx"))
print(summary_tbl, row.names = FALSE)
cat("\nSaved Table_Sensitivity.xlsx\n")
