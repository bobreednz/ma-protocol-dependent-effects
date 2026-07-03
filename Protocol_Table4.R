setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("metafor", "clubSandwich", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(metafor)       # meta-analytic models via rma.mv()
library(clubSandwich)  # CR2 clustered standard errors
library(tidyverse)     # data manipulation
library(writexl)       # export to Excel

# ============================================================
# Load data
# ============================================================

DT <- readRDS("SCData_processed.rds")

# Count total effect sizes and number of studies
n_obs     <- nrow(DT)
n_studies <- length(unique(DT$newid))

# ============================================================
# Impute block-diagonal covariance matrix
# Used by CE and CHE models to account for within-study
# correlation among effect sizes. We assume a constant
# within-study sampling correlation of rho = 0.5.
# ============================================================

rho <- 0.5

V_mat <- vcalc(
  vi      = sez^2,      # sampling variances
  cluster = newid,      # study identifier
  obs     = obs,        # observation identifier within study
  rho     = rho,        # assumed within-study correlation
  data    = DT
)

# ============================================================
# Fit models using REML
# REML produces less biased variance component estimates
# than ML, so we use it for the main reported results.
# ============================================================

# FE: Fixed effects model
# No random component; each effect size is weighted by its
# known sampling variance. Assumes no heterogeneity.
model_FE <- rma.mv(z ~ 1,
                   V      = sez^2,
                   data   = DT,
                   method = "FE")

# RE: Random effects model
# Single variance component at the observation level,
# treating all effect sizes as independently drawn from
# a common distribution.
model_RE <- rma.mv(z ~ 1,
                   V      = sez^2,
                   random = ~ 1 | obs,
                   data   = DT,
                   method = "REML")

# CE: Correlated effects model
# Single between-study variance component, combined with
# the imputed covariance matrix to account for within-study
# correlation among effect sizes.
model_CE <- rma.mv(z ~ 1,
                   V      = V_mat,
                   random = ~ 1 | newid,
                   data   = DT,
                   method = "REML")

# HE: Hierarchical effects model
# Two variance components: between-study (tau) and
# within-study (omega). Effect sizes are treated as
# independent (diagonal V matrix).
model_HE <- rma.mv(z ~ 1,
                   V      = sez^2,
                   random = ~ 1 | newid/obs,
                   data   = DT,
                   method = "REML")

# CHE: Correlated hierarchical effects model
# Two variance components (tau and omega) AND the imputed
# covariance matrix. Most flexible of the five models.
model_CHE <- rma.mv(z ~ 1,
                    V      = V_mat,
                    random = ~ 1 | newid/obs,
                    data   = DT,
                    method = "REML")

# ============================================================
# Refit RE, CE, HE, CHE using ML
# ML log-likelihoods are comparable across models; REML
# log-likelihoods are not. We use these ML versions only
# for likelihood ratio tests.
# ============================================================

model_RE_ML  <- rma.mv(z ~ 1,
                        V      = sez^2,
                        random = ~ 1 | obs,
                        data   = DT,
                        method = "ML")

model_CE_ML  <- rma.mv(z ~ 1,
                        V      = V_mat,
                        random = ~ 1 | newid,
                        data   = DT,
                        method = "ML")

model_HE_ML  <- rma.mv(z ~ 1,
                        V      = sez^2,
                        random = ~ 1 | newid/obs,
                        data   = DT,
                        method = "ML")

model_CHE_ML <- rma.mv(z ~ 1,
                        V      = V_mat,
                        random = ~ 1 | newid/obs,
                        data   = DT,
                        method = "ML")

# ============================================================
# CR2 clustered standard errors
# We use CR2 standard errors throughout, clustered at the
# study level (newid), regardless of which working model
# is used. This makes inference robust to misspecification
# of the variance-covariance structure.
# ============================================================

ct_FE  <- coef_test(model_FE,  cluster = DT$newid, vcov = "CR2")
ct_RE  <- coef_test(model_RE,  cluster = DT$newid, vcov = "CR2")
ct_CE  <- coef_test(model_CE,  cluster = DT$newid, vcov = "CR2")
ct_HE  <- coef_test(model_HE,  cluster = DT$newid, vcov = "CR2")
ct_CHE <- coef_test(model_CHE, cluster = DT$newid, vcov = "CR2")

# ============================================================
# Likelihood ratio tests (using ML models)
#
# We test three nested pairs:
#   (1) RE vs FE: does adding a random effect improve fit?
#   (2) HE vs RE: does a second variance component help?
#   (3) CHE vs CE: does a second variance component help
#       when within-study correlation is already modeled?
#
# CE and HE are not nested with each other (both have one
# variance component but different structures), so no LR
# test is possible between them. Similarly, HE and CHE use
# different covariance matrices (diagonal vs. imputed
# block-diagonal), so no LR test is possible; informal
# evidence comes from the rho profile likelihood (Section 4
# of the protocol), since HE is CHE with rho = 0.
# ============================================================

lrt_RE  <- anova(model_FE,    model_RE_ML)   # RE vs FE
lrt_HE  <- anova(model_RE_ML, model_HE_ML)   # HE vs RE
lrt_CHE <- anova(model_CE_ML, model_CHE_ML)  # CHE vs CE

# ============================================================
# Helper functions for formatting the table
# ============================================================

# Add significance stars based on p-value
stars <- function(p) {
  case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.10 ~ "*",
    TRUE     ~ ""
  )
}

# Format a p-value for the LR test row
fmt_p <- function(p) {
  if (p < 0.0001) "p < .0001"
  else sprintf("p = %.3f", p)
}

# Format a variance component: take sqrt (to get SD) and round
fmt_vc <- function(sigma2) {
  formatC(sqrt(sigma2), digits = 3, format = "f")
}

# Format an estimate with significance stars
fmt_est <- function(ct) {
  paste0(formatC(ct$beta, digits = 3, format = "f"), stars(ct$p_Satt))
}

# Format a standard error in parentheses
fmt_se <- function(ct) {
  paste0("(", formatC(ct$SE, digits = 3, format = "f"), ")")
}

# ============================================================
# Extract LR test p-values
# ============================================================

p_lrt_RE  <- fmt_p(lrt_RE$pval)
p_lrt_HE  <- fmt_p(lrt_HE$pval)
p_lrt_CHE <- fmt_p(lrt_CHE$pval)

# ============================================================
# Assemble the output table
# ============================================================

table_out <- tibble(
  Variable = c("Constant", "(SE)", "Observations", "Studies",
               "tau", "omega", "LR test"),

  # FE: no random effects, so no variance components or LR test
  FE = c(
    fmt_est(ct_FE),
    fmt_se(ct_FE),
    n_obs,
    n_studies,
    "--",
    "--",
    "--"
  ),

  # RE: one variance component (tau); LR test vs FE
  RE = c(
    fmt_est(ct_RE),
    fmt_se(ct_RE),
    n_obs,
    n_studies,
    fmt_vc(model_RE$sigma2),
    "--",
    p_lrt_RE
  ),

  # CE: one between-study variance component (tau) + imputed V;
  # no LR test (not nested with RE or HE)
  CE = c(
    fmt_est(ct_CE),
    fmt_se(ct_CE),
    n_obs,
    n_studies,
    fmt_vc(model_CE$sigma2),
    "--",
    "--"
  ),

  # HE: two variance components (tau, omega); LR test vs RE
  HE = c(
    fmt_est(ct_HE),
    fmt_se(ct_HE),
    n_obs,
    n_studies,
    fmt_vc(model_HE$sigma2[1]),
    fmt_vc(model_HE$sigma2[2]),
    p_lrt_HE
  ),

  # CHE: two variance components (tau, omega) + imputed V;
  # LR test vs CE
  CHE = c(
    fmt_est(ct_CHE),
    fmt_se(ct_CHE),
    n_obs,
    n_studies,
    fmt_vc(model_CHE$sigma2[1]),
    fmt_vc(model_CHE$sigma2[2]),
    p_lrt_CHE
  )
)

# ============================================================
# Export to Excel
# ============================================================

write_xlsx(table_out, "Table4_Protocol.xlsx")

cat("Done. Output saved to Table4_Protocol.xlsx\n")

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
                                                                                                                                