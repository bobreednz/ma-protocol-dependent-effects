# =============================================================================
# Protocol_RhoSensitivity.R
# =============================================================================
# Purpose: Assess the sensitivity of the CHE results to the assumed
# within-study sampling correlation rho.
#
# Throughout the protocol, V_mat is built by vcalc() with an assumed
# within-study correlation rho = 0.5. This value is not identified from the
# data: the true within-study correlations are not observed, and rho enters the
# assumed sampling-covariance matrix rather than the random-effects model. It
# therefore cannot be selected by the likelihood or by any other model-
# comparison criterion, because different values of rho correspond to different
# assumed data rather than to competing models fit to the same data. The
# appropriate check is instead a sensitivity analysis: assume rho = 0.5, then
# confirm that the substantive conclusions do not depend on that choice.
#
# For each candidate rho on a grid:
#   1. Build V_mat with that rho.
#   2. Fit the intercept-only CHE model by REML, which estimates tau^2 and
#      omega^2 at that fixed rho.
#   3. Record the CHE intercept, its CR2 standard error, and tau and omega.
#
# Two figures are produced:
#   Figure_RhoIntercept.png          -- CHE intercept (with 95% CI band) vs rho
#   Figure_RhoVarianceComponents.png -- tau and omega (variance components) vs rho
#
# Both figures include a vertical dashed line at rho = 0.5 (the protocol
# assumption) so the reader can see at a glance whether the intercept and the
# variance components are sensitive to the choice of rho.
#
# Output:
#   Table_RhoSensitivity.xlsx        -- full results for every rho on the grid
#                                        (retained for verification; not
#                                        displayed in the rendered protocol,
#                                        which shows the two figures instead)
#   Figure_RhoIntercept.png
#   Figure_RhoVarianceComponents.png
# =============================================================================

# -- Packages ------------------------------------------------------------------

library(here)           # resolves file paths relative to the project root (.Rproj)
library(metafor)       # rma.mv() and vcalc()
library(clubSandwich)  # coef_test() for CR2 standard errors
library(tidyverse)     # ggplot2 and data manipulation
library(writexl)       # write_xlsx()

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# -- Data ----------------------------------------------------------------------

dat <- readRDS(here("SCData_processed.rds"))

cat(sprintf("Dataset: %d observations from %d studies.\n",
            nrow(dat), length(unique(dat$newid))))

# =============================================================================
# Grid of rho values
# =============================================================================
#
# We try rho from 0 to 0.95 in steps of 0.05, giving 20 candidate values.
# rho = 0 corresponds to the hierarchical effects (HE) model where effect
# sizes within a study are treated as independent after accounting for the
# study-level random effect. rho = 1 is not identified (all within-study
# variation collapses into a single block), so we stop at 0.95.
#
# 20 model fits take roughly 6 minutes in total; no bootstrap is required.

rho_values <- seq(0, 0.95, by = 0.05)
n_rho      <- length(rho_values)

cat(sprintf("\nFitting CHE model at %d values of rho (0 to 0.95).\n", n_rho))

# Pre-allocate results storage.
results <- data.frame(
  rho       = rho_values,
  intercept = NA_real_,   # intercept estimate (Fisher's z)
  se        = NA_real_,   # CR2 standard error of the intercept
  tau       = NA_real_,   # between-study SD (sqrt of sigma2[1])
  omega     = NA_real_    # within-study SD  (sqrt of sigma2[2])
)

# =============================================================================
# Sensitivity loop
# =============================================================================

for (i in seq_len(n_rho)) {

  rho_i <- rho_values[i]

  # Build the imputed block-diagonal covariance matrix for this rho.
  # Within each study block, the off-diagonal entries are rho * sez_j * sez_k,
  # so increasing rho assumes stronger within-study correlation.
  V_i <- vcalc(
    vi      = sez^2,
    cluster = newid,
    obs     = obs,
    rho     = rho_i,
    data    = dat
  )

  # Fit the intercept-only CHE model by REML.
  # random = ~ 1 | newid/obs gives:
  #   sigma2[1] = tau^2  (between-study variance)
  #   sigma2[2] = omega^2 (within-study variance beyond sampling error)
  # We wrap in tryCatch() in case a boundary value causes a convergence issue.
  che_i <- tryCatch(
    rma.mv(
      yi     = z,
      V      = V_i,
      random = ~ 1 | newid/obs,
      data   = dat,
      method = "REML"
    ),
    error   = function(e) {
      cat(sprintf("  ERROR at rho = %.2f: %s\n", rho_i, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(che_i)) {

    # Extract the intercept and its CR2 SE (clustered at the study level).
    ct_i <- coef_test(
      obj     = che_i,
      vcov    = "CR2",
      cluster = dat$newid
    )
    results$intercept[i] <- ct_i[1, "beta"]
    results$se[i]        <- ct_i[1, "SE"]

    # Variance components: sigma2[1] is tau^2, sigma2[2] is omega^2.
    # We report standard deviations (square roots) for interpretability.
    results$tau[i]   <- sqrt(che_i$sigma2[1])
    results$omega[i] <- sqrt(che_i$sigma2[2])

    cat(sprintf("  rho = %.2f  intercept = %.4f  se = %.4f  tau = %.4f  omega = %.4f\n",
                rho_i, results$intercept[i], results$se[i],
                results$tau[i], results$omega[i]))
  }
}

# =============================================================================
# Figure: CHE intercept (with 95% CI band) vs rho (Figure_RhoIntercept.png)
# =============================================================================
#
# The CI band uses intercept +/- 1.96 * SE. This is an approximation because
# the exact critical value uses Satterthwaite degrees of freedom, which vary
# by rho. For a sensitivity display this approximation is adequate.

results_plot <- results %>%
  filter(!is.na(intercept)) %>%
  mutate(
    ci_lo = intercept - 1.96 * se,
    ci_hi = intercept + 1.96 * se
  )

p_int <- ggplot(results_plot, aes(x = rho, y = intercept)) +

  # 95% CI shaded band
  geom_ribbon(
    aes(ymin = ci_lo, ymax = ci_hi),
    fill  = "grey80",
    alpha = 0.6
  ) +

  # Intercept curve
  geom_line(linewidth = 0.9, color = "black") +
  geom_point(size = 2.5, color = "black") +

  # Vertical line at rho = 0.5 (protocol assumption)
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    color      = "grey50",
    linewidth  = 0.7
  ) +
  annotate("text", x = 0.5,
           y = min(results_plot$ci_lo, na.rm = TRUE),
           label = "rho = 0.5\n(protocol)", hjust = -0.1, vjust = 0,
           size = 3, color = "grey40") +

  labs(
    x        = "Within-study correlation (rho)",
    y        = "CHE intercept (Fisher's z)",
    title    = "Sensitivity of the CHE intercept to rho",
    subtitle = "Shaded band = approximate 95% CI (intercept +/- 1.96 * CR2 SE)."
  ) +
  scale_x_continuous(breaks = seq(0, 0.95, by = 0.1)) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9))

ggsave(
  filename = here("Figure_RhoIntercept.png"),
  plot     = p_int,
  width    = 7,
  height   = 5,
  dpi      = 300
)
cat("\nFigure_RhoIntercept.png saved.\n")

# =============================================================================
# Figure: Variance components (tau and omega) vs rho (Figure_RhoVarianceComponents.png)
# =============================================================================
#
# Tau (between-study SD) and omega (within-study SD) are not shown in the
# figure above, so they are plotted together here. As rho increases, more of
# the observed similarity among a study's effect sizes is attributed to
# correlated sampling error rather than to a shared study-level true effect,
# so heterogeneity is mechanically reallocated away from the between-study
# component (tau) and toward the within-study component (omega). Reshape to
# long format so both series share one legend.

results_var <- results %>%
  filter(!is.na(tau)) %>%
  select(rho, tau, omega) %>%
  pivot_longer(
    cols      = c(tau, omega),
    names_to  = "component",
    values_to = "value"
  ) %>%
  mutate(
    component = recode(
      component,
      tau   = "Tau (between-study SD)",
      omega = "Omega (within-study SD)"
    )
  )

p_var <- ggplot(results_var, aes(x = rho, y = value, color = component, linetype = component)) +

  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +

  # Vertical line at rho = 0.5 (protocol assumption)
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    color      = "grey50",
    linewidth  = 0.7
  ) +
  annotate("text", x = 0.5, y = min(results_var$value, na.rm = TRUE),
           label = "rho = 0.5\n(protocol)", hjust = -0.1, vjust = 0,
           size = 3, color = "grey40") +

  scale_color_manual(values = c(
    "Tau (between-study SD)"   = "black",
    "Omega (within-study SD)"  = "#d6604d"
  )) +
  scale_linetype_manual(values = c(
    "Tau (between-study SD)"   = "solid",
    "Omega (within-study SD)"  = "dashed"
  )) +

  labs(
    x        = "Within-study correlation (rho)",
    y        = "Standard deviation (Fisher's z)",
    color    = NULL,
    linetype = NULL,
    title    = "Sensitivity of variance components to rho",
    subtitle = "Between-study (tau) vs. within-study (omega) heterogeneity, CHE model (intercept only)."
  ) +
  scale_x_continuous(breaks = seq(0, 0.95, by = 0.1)) +
  theme_bw(base_size = 11) +
  theme(
    plot.subtitle   = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(
  filename = here("Figure_RhoVarianceComponents.png"),
  plot     = p_var,
  width    = 7,
  height   = 5,
  dpi      = 300
)
cat("Figure_RhoVarianceComponents.png saved.\n")

# =============================================================================
# Results table
# =============================================================================

tbl <- results %>%
  mutate(
    across(c(intercept, se, tau, omega), ~ round(.x, 4))
  ) %>%
  rename(
    `Rho`           = rho,
    `Intercept`     = intercept,
    `SE (CR2)`      = se,
    `Tau`           = tau,
    `Omega`         = omega
  )

write_xlsx(
  list("Rho sensitivity" = tbl),
  path = here("Table_RhoSensitivity.xlsx")
)
cat("Table_RhoSensitivity.xlsx saved.\n")

cat("\nDone. Summary:\n")
cat(sprintf("  Intercept at rho = 0.50:      %.4f (SE = %.4f)\n",
            results$intercept[results$rho == 0.50],
            results$se[results$rho == 0.50]))
cat(sprintf("  Intercept at rho = 0.00:      %.4f\n", results$intercept[results$rho == 0.00]))
cat(sprintf("  Intercept at rho = 0.95:      %.4f\n", results$intercept[results$rho == 0.95]))

# -- Run time ------------------------------------------------------------------
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
