# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("here", "metafor", "tidyverse", "openxlsx")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)       # resolves file paths relative to the project root (.Rproj)
library(metafor)   # rma(), rma.mv(), weights()
library(tidyverse) # dplyr, tibble
library(openxlsx)  # createWorkbook(), setColWidths(), mergeCells()

# ============================================================
# Load data
# ============================================================

DT <- readRDS(here("SCData_processed.rds"))

n_studies <- length(unique(DT$newid))

# Add a unique observation-level ID -- required for the
# three-level HE model, which needs to distinguish
# observations nested within studies.
DT <- DT %>% mutate(obs_id = row_number())

# ============================================================
# FE and RE: fit to all estimates, then sum weights by study
#
# These are the FE and RE models of Table 4 (tbl-table4),
# fitted to all 957 estimates. Each treats every estimate as
# independent. rma() gives the same estimates as the rma.mv()
# calls in Protocol_FiveEstimators.R (RE: random = ~ 1 | obs)
# and its weights() are simpler to extract.
#
# weights() returns each estimate's share of the total weight,
# in percent. Summing within newid gives each study's share.
# Under FE, a study's weight is the sum of the precisions of its
# estimates. Under RE, each estimate's weight is 1/(sez^2 + tau^2),
# so a study's weight still grows with every estimate it reports.
# ============================================================

model_FE <- rma(yi = z, sei = sez, data = DT, method = "FE")
model_RE <- rma(yi = z, sei = sez, data = DT, method = "REML")

study_w <- DT %>%
  mutate(w_FE_obs = as.numeric(weights(model_FE)),
         w_RE_obs = as.numeric(weights(model_RE))) %>%
  group_by(newid) %>%
  summarise(w_FE = sum(w_FE_obs), w_RE = sum(w_RE_obs), .groups = "drop") %>%
  mutate(w_FE = w_FE / sum(w_FE) * 100,   # normalise to exactly 100
         w_RE = w_RE / sum(w_RE) * 100)

# Diagnostic: the RE mean should match Table 4 (0.145)
cat(sprintf("RE mean (all estimates): %.3f\n", as.numeric(coef(model_RE))))

# ============================================================
# HE: Fit three-level model on full data, extract GLS weights
#
# The HE (Hierarchical Effects) estimator uses a three-level
# random effects structure:
#   Level 1: sampling error         (obs_id within newid)
#   Level 2: within-study variance  (sigma2[2])
#   Level 3: between-study variance (sigma2[1])
#
# Two key implementation points:
#
#   (1) We fit on the full dataset (not study-aggregated).
#       Effects within the same study are correlated through
#       the shared between-study variance component, so the
#       correct GLS weights require inverting the full
#       block-diagonal covariance matrix -- not treating
#       observations as independent.
#
#   (2) We extract weights with weights(model_HE, type = "rowsum"),
#       which returns the row sums of M^{-1} -- the correct GLS
#       weights for an intercept-only model. The default
#       type = "diagonal" returns only the diagonal elements of
#       M^{-1}, which ignores the negative off-diagonal terms
#       within studies and over-weights large studies. We sum
#       the row-sum weights within each study to get study-level
#       weights.
# ============================================================

model_HE <- rma.mv(
  yi     = z,
  V      = sez^2,          # sampling variance at observation level
  random = ~ 1 | newid / obs_id,   # between-study / within-study
  data   = DT,
  method = "REML"
)

# For an intercept-only rma.mv model, the correct GLS weights are
# the ROW SUMS of M^{-1} (the inverse marginal variance-covariance
# matrix), not the diagonal elements. The row sums account for the
# negative off-diagonal terms within studies, which reduce the
# effective weight of observations from the same study.
# type = "rowsum" returns these correct weights as percentages
# summing to 100. The default (type = "diagonal") gives only the
# diagonal elements and over-weights studies with many observations.
w_obs_HE <- weights(model_HE, type = "rowsum")

# Diagnostic: verify weights sum to 100 and reproduce the pooled estimate.
stopifnot(abs(sum(w_obs_HE) - 100) < 1e-4)
mu_check <- sum((w_obs_HE / 100) * DT$z)
stopifnot(isTRUE(all.equal(as.numeric(coef(model_HE)), mu_check, tolerance = 1e-4)))

# Sum observation-level weights within studies (preserving newid type).
study_w_HE <- DT %>%
  mutate(w_HE_obs = as.numeric(w_obs_HE)) %>%
  group_by(newid) %>%
  summarise(w_HE = sum(w_HE_obs), .groups = "drop") %>%
  mutate(w_HE = w_HE / sum(w_HE) * 100)  # normalise to exactly 100

study_w <- study_w %>% left_join(study_w_HE, by = "newid")

# ============================================================
# I-squared and tau for RE
#
# RE is fitted to all estimates without study structure, so its
# tau is the total heterogeneity across estimates, as in Table 4.
# We report tau (the SD, sqrt of tau-squared) rather than
# tau-squared since tau is on the same scale as the effect sizes.
# No CIs are reported for any heterogeneity statistics.
# ============================================================

tau2_est <- model_RE$tau2
tau_RE   <- sqrt(tau2_est)   # total SD of true effects across estimates
i2_est   <- model_RE$I2

# ============================================================
# I-squared, tau (study), and tau (estimate) for HE
#
# Notation follows the multilevel meta-analysis literature:
#   tau (study)    = sqrt(sigma2[1]) = between-study SD
#   tau (estimate) = sqrt(sigma2[2]) = within-study SD (between-estimate)
#
# I-squared components use the "typical" sampling variance
# (Higgins and Thompson, 2002; Cheung, 2014; Konstantopoulos, 2011).
# Profile-likelihood CIs for rma.mv are very slow to compute,
# so only point estimates are reported for HE.
# ============================================================

tau2_between <- model_HE$sigma2[1]  # between-study variance
tau2_within  <- model_HE$sigma2[2]  # within-study variance

tau_between  <- sqrt(tau2_between)  # between-study SD (tau)
omega_within <- sqrt(tau2_within)   # within-study SD (omega)
tau_total_HE <- sqrt(tau2_between + tau2_within)  # total SD, comparable to RE

# Typical sampling variance (Higgins and Thompson, 2002)
vi_vec <- DT$sez^2
w_vec  <- 1 / vi_vec
v_bar  <- ((length(vi_vec) - 1) * sum(w_vec)) / (sum(w_vec)^2 - sum(w_vec^2))

total_var   <- v_bar + tau2_within + tau2_between
I2_between  <- tau2_between / total_var * 100
I2_within   <- tau2_within  / total_var * 100
I2_total_HE <- (tau2_between + tau2_within) / total_var * 100

# ============================================================
# Helper functions
# ============================================================

# Format a weight as a percentage with one decimal place
fmt_pct <- function(x) paste0(formatC(x, digits = 1, format = "f"), "%")

# Percentile shorthand
pct <- function(x, p) quantile(x, p / 100, names = FALSE)

# ============================================================
# Compute top-3 and top-10 weights and identify those studies
# ============================================================

study_w_sorted_FE <- study_w[order(study_w$w_FE, decreasing = TRUE), ]
study_w_sorted_RE <- study_w[order(study_w$w_RE, decreasing = TRUE), ]
study_w_sorted_HE <- study_w[order(study_w$w_HE, decreasing = TRUE), ]

top3_FE   <- sum(study_w_sorted_FE$w_FE[1:3])
top10_FE  <- sum(study_w_sorted_FE$w_FE[1:10])
top3_RE   <- sum(study_w_sorted_RE$w_RE[1:3])
top10_RE  <- sum(study_w_sorted_RE$w_RE[1:10])
top3_HE   <- sum(study_w_sorted_HE$w_HE[1:3])
top10_HE  <- sum(study_w_sorted_HE$w_HE[1:10])

top3_ids  <- study_w_sorted_FE$newid[1:3]
top10_ids <- study_w_sorted_FE$newid[1:10]

cat(sprintf("Top 3 studies by FE weight (id): %s\n",
            paste(top3_ids, collapse = ", ")))
cat(sprintf("Top 10 studies by FE weight (id): %s\n",
            paste(top10_ids, collapse = ", ")))

# ============================================================
# Assemble the table
# ============================================================

table3 <- tibble(
  # tau (study) = between-study SD; tau (estimate) = within-study SD
  # (HE only); tau (total) = square root of the sum of the variances.
  # I-squared comes last: it depends on the precision of the
  # estimates, so it is not a measure of how much the true effects vary.
  Statistic = c("Mean", "Median", "5%", "10%", "90%", "95%",
                "Minimum", "Maximum", "Top 3", "Top 10",
                "tau (study)", "tau (estimate)", "tau (total)", "I-squared", "Studies"),

  FE = c(
    fmt_pct(mean(study_w$w_FE)),
    fmt_pct(median(study_w$w_FE)),
    fmt_pct(pct(study_w$w_FE,  5)),
    fmt_pct(pct(study_w$w_FE, 10)),
    fmt_pct(pct(study_w$w_FE, 90)),
    fmt_pct(pct(study_w$w_FE, 95)),
    fmt_pct(min(study_w$w_FE)),
    fmt_pct(max(study_w$w_FE)),
    fmt_pct(top3_FE),
    fmt_pct(top10_FE),
    "--",   # tau (study)
    "--",   # tau (estimate)
    "--",   # tau (total)
    "--",   # I-squared
    as.character(n_studies)
  ),

  RE = c(
    fmt_pct(mean(study_w$w_RE)),
    fmt_pct(median(study_w$w_RE)),
    fmt_pct(pct(study_w$w_RE,  5)),
    fmt_pct(pct(study_w$w_RE, 10)),
    fmt_pct(pct(study_w$w_RE, 90)),
    fmt_pct(pct(study_w$w_RE, 95)),
    fmt_pct(min(study_w$w_RE)),
    fmt_pct(max(study_w$w_RE)),
    fmt_pct(top3_RE),
    fmt_pct(top10_RE),
    # RE is fitted to all estimates without study structure, so its
    # single tau is the total heterogeneity across estimates
    "--",                                        # tau (study)
    "--",                                        # tau (estimate)
    formatC(tau_RE, digits = 4, format = "f"),   # tau (total), no CI
    paste0(formatC(i2_est, digits = 1, format = "f"), "%"),
    as.character(n_studies)
  ),

  HE = c(
    fmt_pct(mean(study_w$w_HE)),
    fmt_pct(median(study_w$w_HE)),
    fmt_pct(pct(study_w$w_HE,  5)),
    fmt_pct(pct(study_w$w_HE, 10)),
    fmt_pct(pct(study_w$w_HE, 90)),
    fmt_pct(pct(study_w$w_HE, 95)),
    fmt_pct(min(study_w$w_HE)),
    fmt_pct(max(study_w$w_HE)),
    fmt_pct(top3_HE),
    fmt_pct(top10_HE),
    formatC(tau_between,  digits = 4, format = "f"),  # tau (study): between-study SD
    formatC(omega_within, digits = 4, format = "f"),  # tau (estimate): within-study SD
    formatC(tau_total_HE, digits = 4, format = "f"),  # tau (total)
    paste0(formatC(I2_total_HE, digits = 1, format = "f"), "%"),
    as.character(n_studies)
  )
)

# ============================================================
# Export to Excel
#
# Footnotes are NOT included in the Excel file. read_excel()
# converts empty cells to NA, so mixed-width footnote rows
# cannot be stored cleanly in a rectangular table. Instead,
# footnotes are printed to the console here (to be copied
# into the Quarto document as plain text below the table).
#
# openxlsx is used for column width and wrap-text control.
# ============================================================

cat("\nTable 3 footnotes (add to Quarto document below the table):\n")
cat(sprintf("a  Studies with the three largest FE weights (id): %s.\n",
            paste(top3_ids, collapse = ", ")))
cat(sprintf("b  Studies with the ten largest FE weights (id): %s.\n",
            paste(top10_ids, collapse = ", ")))

wb <- createWorkbook()
addWorksheet(wb, "Table3")

writeData(wb, "Table3", table3, startRow = 1, startCol = 1,
          headerStyle = createStyle(textDecoration = "bold"))

# Column widths: Statistic | FE | RE | HE
setColWidths(wb, "Table3", cols = 1:4, widths = c(16, 10, 12, 42))

# Wrap text and top-align all cells
wrap_style <- createStyle(wrapText = TRUE, valign = "top")
addStyle(wb, "Table3", style = wrap_style,
         rows = 2:(nrow(table3) + 1), cols = 1:4,
         gridExpand = TRUE)

saveWorkbook(wb, here("Table_StudyWeights.xlsx"), overwrite = TRUE)
cat("Table 3 saved to Table_StudyWeights.xlsx\n")

# ── Run time ──────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
