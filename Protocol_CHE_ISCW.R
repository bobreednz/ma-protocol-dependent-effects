# =============================================================================
# Protocol_CHE_ISCW.R
# =============================================================================
# Purpose: robustness check on the weights used by the PET-PEESE correction
# (Section 5.2), using the inverse sampling covariance weights (ISCW) of
# Chen and Pustejovsky (2026).
#
# The robustness table has two columns, both estimated without moderators:
#   1. CHE-PET-PEESE: the main PET-PEESE estimate, repeated for comparison.
#      Intercept plus bias term, conventional CHE weights (the inverse of the
#      total variance: sampling error plus the CHE random effects).
#   2. CHE-ISCW-PET-PEESE: intercept plus bias term, with ISCW weights. ISCW
#      weights each study by the inverse of its sampling covariance matrix
#      alone, which gives more weight to large studies. It uses the SAME bias
#      term (sez or sez^2) that the PET-PEESE rule selected for column 1, so
#      the two columns differ only in their weights.
#
# PET-PEESE rule (applied to the conventional CHE model only): report PEESE
# (bias term sez^2) if the PET intercept is positive and its two-sided p-value
# is below 0.10 (a one-sided test at 5 percent); otherwise report PET (bias
# term sez). This is the rule used in the code of Chen and Pustejovsky (2026).
#
# The table also reports the Satterthwaite degrees of freedom of the mean
# effect. ISCW weights can concentrate on a few large studies, which leaves
# few degrees of freedom; Tipton (2015) advises against interpreting RVE
# estimates with fewer than 4.
#
# Output:
#   Table_CHE_ISCW.xlsx  "Display": the two-column robustness table
#                        "Values":  numeric values used by
#                                   Protocol_Diagnostics.R
# Run time: about a minute.
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

# ---- Packages ----------------------------------------------------------------
packages <- c("here", "metafor", "clubSandwich", "writexl")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)          # resolves file paths relative to the project root (.Rproj)
library(metafor)       # rma.mv() and vcalc()
library(clubSandwich)  # CR2 cluster-robust standard errors (used by robust())
library(writexl)       # write_xlsx()

# ---- Data --------------------------------------------------------------------
dat <- readRDS(here("SCData_processed.rds"))
cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# ---- Sampling covariance matrix and ISCW weights ---------------------------
# Same imputed covariance matrix as every CHE model in the protocol (rho = 0.5).
V_mat <- vcalc(vi = sez^2, cluster = newid, obs = obs, rho = 0.5, data = dat)

# The ISCW weight matrix is simply the inverse of the sampling covariance matrix.
W_iscw <- solve(V_mat)

# ---- Helper: fit one CHE model and return its estimates ---------------------
# rhs:  right-hand side of the formula (the bias term)
# iscw: TRUE uses ISCW weights; FALSE uses the conventional CHE weights
fit_che <- function(rhs, iscw) {
  args <- list(yi = dat$z, V = V_mat, mods = as.formula(paste("~", rhs)),
               random = ~ 1 | newid/obs, data = dat, method = "REML")
  if (iscw) args$W <- W_iscw
  m   <- do.call(rma.mv, args)
  rob <- robust(m, cluster = dat$newid, clubSandwich = TRUE)   # CR2 inference
  list(est = as.numeric(rob$beta[1]), se = rob$se[1], p = rob$pval[1],
       lower = rob$ci.lb[1], upper = rob$ci.ub[1],
       df = rob$ddf[1],                     # Satterthwaite df of the intercept
       bias = as.numeric(rob$beta[2]), bias_se = rob$se[2],
       tau = sqrt(m$sigma2[1]), omega = sqrt(m$sigma2[2]))
}

# The two bias terms: PET uses the standard error, PEESE its square.
bias_term <- c(PET = "sez", PEESE = "I(sez^2)")

# Apply the PET-PEESE rule to the conventional CHE model: fit PET and PEESE,
# keep PEESE only if the PET intercept is positive with a two-sided p-value
# below 0.10.
cat("\nFitting CHE-PET-PEESE (conventional weights)...\n")
pet   <- fit_che(bias_term["PET"],   iscw = FALSE)
peese <- fit_che(bias_term["PEESE"], iscw = FALSE)
use_peese <- pet$est > 0 && pet$p < 0.10
che_pp <- if (use_peese) peese else pet
che_pp$correction <- if (use_peese) "PEESE" else "PET"

# Fit the ISCW version with the bias term already selected above. The rule is
# NOT reapplied, so the only difference from the conventional model is the
# weights.
cat(sprintf("Fitting CHE-ISCW-PET-PEESE (bias term from CHE-PET-PEESE: %s)...\n",
            che_pp$correction))
iscw_pp <- fit_che(bias_term[che_pp$correction], iscw = TRUE)
iscw_pp$correction <- che_pp$correction

# ---- Build the tables ----------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
ci <- function(lo, hi) sprintf("[%.3f, %.3f]", lo, hi)
column <- function(r) c(
  f3(r$est), paste0("(", f3(r$se), ")"), ci(r$lower, r$upper),
  sprintf("%.1f", r$df),
  f3(r$bias), paste0("(", f3(r$bias_se), ")"),
  f3(r$tau), f3(r$omega), nrow(dat), length(unique(dat$newid)))
row_labels <- c("Mean effect", "(SE)", "95% CI", "Degrees of freedom",
                "Bias coefficient", "(SE)", "tau", "omega", "Observations", "Studies")

display <- data.frame(
  Parameter = row_labels,
  `CHE-PET-PEESE` = column(che_pp),
  `CHE-ISCW-PET-PEESE` = column(iscw_pp),
  check.names = FALSE, stringsAsFactors = FALSE)

values <- data.frame(
  Quantity   = c("CHE-PET-PEESE", "CHE-ISCW-PET-PEESE"),
  Correction = c(che_pp$correction, iscw_pp$correction),
  Estimate   = c(che_pp$est, iscw_pp$est),
  Lower      = c(che_pp$lower, iscw_pp$lower),
  Upper      = c(che_pp$upper, iscw_pp$upper),
  df         = c(che_pp$df, iscw_pp$df),
  stringsAsFactors = FALSE)

write_xlsx(list(Display = display, Values = values),
           path = here("Table_CHE_ISCW.xlsx"))
print(display, row.names = FALSE)
cat("\nSaved Table_CHE_ISCW.xlsx\n")

elapsed <- proc.time() - start_time
cat(sprintf("Total run time: %.1f seconds.\n", elapsed["elapsed"]))
