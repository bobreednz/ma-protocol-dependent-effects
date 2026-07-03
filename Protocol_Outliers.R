# =============================================================================
# Protocol_Outliers.R
# =============================================================================
# Purpose: Outlier and influence analysis for the appendix.
#
# This script has two parts:
#
#   Part 1 -- Observation-level screening using influence() on an RE model.
#             influence() is not available for rma.mv() objects, so we fit
#             a standard rma() RE model as a screening device. This ignores
#             within-study dependence but gives a useful first pass at
#             identifying extreme individual effect sizes.
#
#   Part 2 -- Study-level jackknife on the CHE model. We drop one study at
#             a time, refit the full CHE model, and record how the intercept
#             changes. This is the primary robustness check because it (a)
#             uses the correct estimator and (b) operates at the study level,
#             which is the natural unit of dependence in a meta-analysis with
#             clustered effect sizes.
#
# Outputs:
#   Figure_Influence_RE.png     -- influence diagnostics for the RE model
#   Figure_Jackknife_CHE.png    -- leave-one-study-out intercept estimates
#   Table_Jackknife_CHE.xlsx    -- all studies sorted by absolute change
#                                  in the leave-one-out intercept
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(metafor)       # rma() and rma.mv() for meta-analytic models
library(clubSandwich)  # coef_test() for CR2 clustered standard errors
library(tidyverse)     # data manipulation and ggplot2 for figures
library(writexl)       # write_xlsx() to export Excel output

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ── Data ──────────────────────────────────────────────────────────────────────

# Load the processed dataset produced by Protocol_DataCleaning.R.
# Never use the raw SCData (20240508).dta directly.
dat <- readRDS("SCData_processed.rds")

cat(sprintf("Dataset loaded: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# =============================================================================
# Part 1: Observation-level influence diagnostics (RE model)
# =============================================================================
#
# We fit a standard univariate RE model using rma() -- not rma.mv() -- because
# influence() is only defined for rma() objects. Treating all 957 observations
# as independent understates uncertainty, but that does not matter here: we are
# using this model as a screening tool to flag individual effect sizes that are
# extreme relative to the pooled estimate, not as a primary estimator.

re_screen <- rma(
  yi     = z,       # Fisher's z effect size
  sei    = sez,     # standard error of z
  data   = dat,
  method = "REML"
)

# Compute the full suite of influence diagnostics. For each observation the
# output includes: externally studentized residuals (rstudent), DFFITS,
# Cook's distance, covariance ratio, leave-one-out tau^2 and QE, hat value,
# and weight.
inf <- influence(re_screen)

# Save the influence plot. The default plot() method for influence objects
# produces one panel per diagnostic. With 957 observations the plot is dense,
# so we use a tall canvas and small point size.
png(
  filename = "Figure_Influence_RE.png",
  width    = 10,
  height   = 14,
  units    = "in",
  res      = 300
)
plot(
  inf,
  cex = 0.35,   # small points to avoid overplotting at n = 957
  las = 1       # horizontal tick labels on all axes
)
dev.off()

cat("Figure_Influence_RE.png saved.\n")

# Print the observations flagged as extreme by the studentized residual
# criterion (|rstudent| > 3). This is a conventional threshold.
# inf is an "infl.rma.uni" object whose diagnostics live in inf$inf.
# We extract that inner data frame directly; as.data.frame() on the outer
# object fails because the class is not a plain data frame.
inf_df <- as.data.frame(inf$inf)
extreme_obs <- inf_df[abs(inf_df$rstudent) > 3, ]

cat(sprintf("\nObservations with |studentized residual| > 3: %d\n",
            nrow(extreme_obs)))
if (nrow(extreme_obs) > 0) print(extreme_obs)

# =============================================================================
# Part 2: Study-level jackknife on the CHE model
# =============================================================================

# ── 2a. Full-sample CHE model ─────────────────────────────────────────────────

# Build the imputed block-diagonal covariance matrix. Within each study block,
# off-diagonal entries equal rho * sez_i * sez_j (i.e., the two sampling SEs
# times the assumed within-study correlation). The assumed correlation rho = 0.5
# is the same value used throughout the protocol. metafor::vcalc() replaced
# the deprecated clubSandwich::impute_covariance_matrix().
V_full <- vcalc(
  vi      = sez^2,       # sampling variances (SE squared)
  cluster = newid,       # study identifier for blocking
  obs     = obs,         # observation identifier within study
  rho     = 0.5,         # assumed within-study correlation
  data    = dat
)

# Fit the full CHE (correlated-and-hierarchical effects) model.
# The random effects specification ~ 1 | newid/obs gives two variance
# components:
#   sigma2[1] = tau^2  (between-study heterogeneity)
#   sigma2[2] = omega^2 (within-study heterogeneity beyond sampling variance)
che_full <- rma.mv(
  yi     = z,
  V      = V_full,
  random = ~ 1 | newid/obs,
  data   = dat,
  method = "REML"
)

# Extract the intercept and its SE using CR2 clustered standard errors,
# clustered at the study level (newid). CR2 SEs are robust to misspecification
# of the variance structure.
ct_full <- coef_test(
  obj     = che_full,
  vcov    = "CR2",
  cluster = dat$newid
)

# The intercept is always the first row of coef_test() output.
intercept_full <- ct_full[1, "beta"]
se_full        <- ct_full[1, "SE"]

cat(sprintf("\nFull-sample CHE intercept: %.4f (SE = %.4f)\n",
            intercept_full, se_full))

# ── 2b. Jackknife loop ────────────────────────────────────────────────────────

study_ids <- sort(unique(dat$newid))
n_studies <- length(study_ids)

cat(sprintf("\nRunning leave-one-study-out jackknife over %d studies.\n",
            n_studies))
cat("This takes approximately 2 hours (83 CHE model fits); consider running overnight.\n\n")

# Pre-allocate a results data frame with one row per study.
jk <- data.frame(
  study     = study_ids,
  n_obs     = NA_integer_,   # number of effect sizes from the dropped study
  intercept = NA_real_,      # leave-one-out intercept
  se        = NA_real_,      # leave-one-out SE (CR2)
  tau       = NA_real_,      # leave-one-out tau (between-study SD)
  omega     = NA_real_       # leave-one-out omega (within-study SD)
)

for (i in seq_along(study_ids)) {

  sid <- study_ids[i]

  # Count how many effect sizes this study contributes to the full sample.
  jk$n_obs[i] <- sum(dat$newid == sid)

  # Drop study sid and refit on the remaining 82 studies.
  dat_loo <- dat[dat$newid != sid, ]

  # Rebuild the covariance matrix for the reduced dataset. We cannot reuse
  # V_full because the matrix depends on which studies are present; the
  # dropped study's rows must be removed.
  V_loo <- vcalc(
    vi      = sez^2,
    cluster = newid,
    obs     = obs,
    rho     = 0.5,
    data    = dat_loo
  )

  # Refit CHE. tryCatch() ensures that a convergence failure for any single
  # study does not crash the loop -- the row for that study is left as NA.
  che_loo <- tryCatch(
    rma.mv(
      yi     = z,
      V      = V_loo,
      random = ~ 1 | newid/obs,
      data   = dat_loo,
      method = "REML"
    ),
    error = function(e) {
      cat(sprintf("  ERROR for study %d: %s\n", sid, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(che_loo)) {
    ct_loo <- coef_test(
      obj     = che_loo,
      vcov    = "CR2",
      cluster = dat_loo$newid
    )
    jk$intercept[i] <- ct_loo[1, "beta"]
    jk$se[i]        <- ct_loo[1, "SE"]
    # sigma2[1] is tau^2 (between-study); sigma2[2] is omega^2 (within-study).
    # We report SDs (square roots) to match the protocol convention.
    jk$tau[i]   <- sqrt(che_loo$sigma2[1])
    jk$omega[i] <- sqrt(che_loo$sigma2[2])
  }

  # Print progress every 10 studies.
  if (i %% 10 == 0) {
    cat(sprintf("  Completed %d of %d studies\n", i, n_studies))
  }
}

# The "change" column is the difference between the leave-one-out intercept
# and the full-sample intercept. A large positive value means that dropping
# that study raises the overall estimate; a large negative value means it
# lowers it.
jk$change <- jk$intercept - intercept_full

cat("\nJackknife complete.\n")
cat(sprintf("Largest upward shift:   +%.4f (study %d)\n",
            max(jk$change, na.rm = TRUE),
            jk$study[which.max(jk$change)]))
cat(sprintf("Largest downward shift: %.4f (study %d)\n",
            min(jk$change, na.rm = TRUE),
            jk$study[which.min(jk$change)]))

# ── 2c. Jackknife figure ──────────────────────────────────────────────────────

# Build the plot data. We sort by leave-one-out intercept (ascending) so
# that the forest-plot-style display has a clear ordering.
jk_plot <- jk %>%
  filter(!is.na(intercept)) %>%
  arrange(intercept) %>%
  mutate(
    # Create an ordered factor so ggplot respects our sort order.
    study_label = factor(
      paste0("Study ", study),
      levels = paste0("Study ", study)
    ),
    # Approximate 95% CIs using the normal distribution. These are only for
    # visual guidance; formal inference uses Satterthwaite df elsewhere.
    lower = intercept - 1.96 * se,
    upper = intercept + 1.96 * se,
    # Color points by direction: does dropping this study raise or lower
    # the intercept relative to the full-sample estimate?
    direction = ifelse(change >= 0, "Rises", "Falls")
  )

p_jk <- ggplot(jk_plot, aes(x = intercept, y = study_label)) +

  # 95% CI horizontal bars
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",   # horizontal bars (the modern replacement for geom_errorbarh)
    width       = 0.4,
    color       = "grey65",
    linewidth   = 0.35
  ) +

  # Leave-one-out point estimates, colored by direction of change
  geom_point(
    aes(color = direction),
    size = 1.6
  ) +

  # Vertical dashed line at the full-sample estimate
  geom_vline(
    xintercept = intercept_full,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.55
  ) +

  scale_color_manual(
    values = c("Rises" = "#2166ac",   # blue  = estimate rises when study dropped
               "Falls" = "#d6604d"),  # red   = estimate falls when study dropped
    name   = "Effect of dropping study on intercept"
  ) +

  labs(
    x        = "Leave-one-out intercept (Fisher's z)",
    y        = NULL,
    title    = "Study-level jackknife: CHE intercept",
    subtitle = paste0(
      "Dashed line = full-sample estimate (z = ",
      round(intercept_full, 3),
      "). Points sorted by LOO intercept. Bars show approximate 95% CIs."
    )
  ) +

  theme_bw(base_size = 9) +
  theme(
    legend.position    = "bottom",
    axis.text.y        = element_text(size = 6.5),
    panel.grid.major.y = element_blank(),
    plot.subtitle      = element_text(size = 7.5)
  )

ggsave(
  filename = "Figure_Jackknife_CHE.png",
  plot     = p_jk,
  width    = 7,
  height   = 11,
  dpi      = 300
)

cat("Figure_Jackknife_CHE.png saved.\n")

# ── 2d. Companion table ───────────────────────────────────────────────────────

# All 83 studies are included in the table, sorted by absolute change in the
# intercept so the most influential studies appear at the top.
jk_table <- jk %>%
  arrange(desc(abs(change))) %>%
  transmute(
    `Study`            = study,
    `N (obs)`          = n_obs,
    `Full intercept`   = round(intercept_full, 4),
    `LOO intercept`    = round(intercept, 4),
    `LOO SE`           = round(se, 4),
    `Change`           = round(change, 4),
    `LOO tau`          = round(tau, 4),
    `LOO omega`        = round(omega, 4)
  )

write_xlsx(
  list("Jackknife" = jk_table),
  path = "Table_Jackknife_CHE.xlsx"
)

cat("Table_Jackknife_CHE.xlsx saved.\n")

cat("\nAll outputs written:\n")
cat("  Figure_Influence_RE.png\n")
cat("  Figure_Jackknife_CHE.png\n")
cat("  Table_Jackknife_CHE.xlsx\n")

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - sta