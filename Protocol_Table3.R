setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("metafor", "tidyverse", "openxlsx")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(metafor)   # rma(), rma.mv(), weights()
library(tidyverse) # dplyr, tibble
library(openxlsx)  # createWorkbook(), setColWidths(), mergeCells()

# ============================================================
# Load data
# ============================================================

DT <- readRDS("SCData_processed.rds")

n_studies <- length(unique(DT$newid))

# Add a unique observation-level ID -- required for the
# three-level HE model, which needs to distinguish
# observations nested within studies.
DT <- DT %>% mutate(obs_id = row_number())

# ============================================================
# FE and RE: Aggregate to study level, then fit rma()
#
# We use rma() rather than rma.mv() here because:
#   (1) The weights() function works cleanly with rma().
#   (2) The purpose of Table 3 is to illustrate the
#       weighting properties of FE and RE as simple
#       benchmarks before introducing the multi-level
#       models. Using rma() is consistent with the protocol.
# Note: rma() treats each effect size as independent, which
# ignores within-study dependence. This is intentional here
# since the table is about weighting, not inference.
#
# Aggregation:
#   Study-level estimate: precision-weighted mean of z within
#     each study (weights = 1/sez^2)
#   Study-level variance: 1 / sum(1/sez^2) within each study
# ============================================================

study_agg <- DT %>%
  group_by(newid) %>%
  summarise(
    z_study  = sum(z / sez^2) / sum(1 / sez^2),  # precision-weighted mean
    v_study  = 1 / sum(1 / sez^2),                # study-level sampling variance
    .groups  = "drop"
  ) %>%
  mutate(se_study = sqrt(v_study))

# Fit FE and RE at the study level
model_FE <- rma(yi = z_study, sei = se_study, data = study_agg, method = "FE")
model_RE <- rma(yi = z_study, sei = se_study, data = study_agg, method = "REML")

# Study-level weights from FE and RE (already one row per study)
study_w <- study_agg %>%
  mutate(w_FE = weights(model_FE),
         w_RE = weights(model_RE)) %>%
  select(newid, w_FE, w_RE)

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
# We report tau (the SD, sqrt of tau-squared) rather than
# tau-squared since tau is on the same scale as the effect sizes.
# No CIs are reported for any heterogeneity statistics.
# ============================================================

tau2_est <- model_RE$tau2
tau_RE   <- sqrt(tau2_est)   # between-study SD
i2_est   <- model_RE$I2

# ============================================================
# I-squared, tau, and omega for HE
#
# Notation follows the multilevel meta-analysis literature:
#   tau   = sqrt(sigma2[1]) = between-study SD
#   omega = sqrt(sigma2[2]) = within-study SD (between-estimate)
#
# I-squared components use the harmonic mean of the sampling
# variances (Cheung, 2014; Konstantopoulos, 2011).
# Profile-likelihood CIs for rma.mv are very slow to compute,
# so only point estimates are reported for HE.
# ============================================================

tau2_between <- model_HE$sigma2[1]  # between-study variance
tau2_within  <- model_HE$sigma2[2]  # within-study variance

tau_between  <- sqrt(tau2_between)  # between-study SD (tau)
omega_within <- sqrt(tau2_within)   # within-study SD (omega)

# Harmonic mean of the sampling variances
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
  # tau = between-study SD; omega = within-study SD (HE only)
  Statistic = c("Mean", "Median", "5%", "10%", "90%", "95%",
                "Minimum", "Maximum", "Top 3", "Top 10",
                "I-squared", "tau", "omega", "Studies"),

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
    "--",   # I-squared
    "--",   # tau
    "--",   # omega
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
    paste0(formatC(i2_est, digits = 1, format = "f"), "%"),
    formatC(tau_RE, digits = 4, format = "f"),   # tau (no CI)
    "--",                                        # omega (RE has no within-study component)
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
    paste0(formatC(I2_total_HE, digits = 1, format = "f"), "%"),
    formatC(tau_between,  digits = 4, format = "f"),  # tau  (between-study SD)
    formatC(omega_within, digits = 4, format = "f"),  # omega (within-study SD)
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

saveWorkbook(wb, "Table3_Protocol.xlsx", overwrite = TRUE)
cat("Table 3 saved to Table3_Protocol.xlsx\n")

# ── Run time ──────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
