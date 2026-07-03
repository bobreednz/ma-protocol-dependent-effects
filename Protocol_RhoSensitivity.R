# =============================================================================
# Protocol_RhoSensitivity.R
# =============================================================================
# Purpose: Choose the within-study correlation rho for the imputed
# block-diagonal covariance matrix.
#
# Throughout the protocol, V_mat is built by vcalc() with an assumed
# within-study correlation rho = 0.5. This assumption is not
# identified from the data in the usual sense -- the true within-study
# correlations are not observed -- but we can still ask: for what value of rho
# does the CHE model fit the data best?
#
# The approach is a profile REML likelihood over rho. For each candidate rho:
#   1. Build V_mat with that rho.
#   2. Fit the intercept-only CHE model by REML, which estimates tau^2 and
#      omega^2 at that fixed rho.
#   3. Record the REML log-likelihood.
#
# Because the fixed effects structure is the same in every model (intercept
# only), the REML log-likelihoods are directly comparable across rho values.
# The rho that maximises the REML log-likelihood is the best-fitting choice
# given the data.
#
# Two figures are produced:
#   Figure_RhoLikelihood.png  -- REML log-likelihood vs rho
#   Figure_RhoIntercept.png   -- CHE intercept (with 95% CI band) vs rho
#
# Both figures include a vertical dashed line at rho = 0.5 (the protocol
# assumption) and a vertical solid line at the rho that maximises the
# log-likelihood, so the reader can see at a glance whether the assumed value
# is close to optimal and whether the intercept estimate is sensitive to the
# choice.
#
# Output:
#   Table_RhoSensitivity.xlsx   -- full results for every rho on the grid
#   Figure_RhoLikelihood.png
#   Figure_RhoIntercept.png
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

library(metafor)       # rma.mv() and vcalc()
library(clubSandwich)  # coef_test() for CR2 standard errors
library(tidyverse)     # ggplot2 and data manipulation
library(writexl)       # write_xlsx()

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ── Data ──────────────────────────────────────────────────────────────────────

dat <- readRDS("SCData_processed.rds")

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
  logLik    = NA_real_,   # REML log-likelihood at each rho
  intercept = NA_real_,   # intercept estimate (Fisher's z)
  se        = NA_real_,   # CR2 standard error of the intercept
  tau       = NA_real_,   # between-study SD (sqrt of sigma2[1])
  omega     = NA_real_    # within-study SD  (sqrt of sigma2[2])
)

# =============================================================================
# Profile likelihood loop
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

    # Extract the REML log-likelihood. In metafor, rma.mv objects do not have
    # a $logLik slot; the log-likelihood is retrieved via the logLik() function.
    results$logLik[i] <- as.numeric(logLik(che_i))

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

    cat(sprintf("  rho = %.2f  logLik = %8.3f  intercept = %.4f  tau = %.4f  omega = %.4f\n",
                rho_i, results$logLik[i], results$intercept[i],
                results$tau[i], results$omega[i]))
  }
}

# =============================================================================
# Identify the best rho
# =============================================================================

best_idx <- which.max(results$logLik)
best_rho <- results$rho[best_idx]

cat(sprintf("\nBest rho (maximum REML log-likelihood): %.2f\n", best_rho))
cat(sprintf("  logLik at best rho:   %.4f\n", results$logLik[best_idx]))
cat(sprintf("  logLik at rho = 0.50: %.4f\n",
            results$logLik[results$rho == 0.50]))
cat(sprintf("  Difference:           %.4f\n",
            results$logLik[best_idx] - results$logLik[results$rho == 0.50]))

# =============================================================================
# Figure: REML log-likelihood vs rho (Figure_RhoLikelihood.png)
# =============================================================================

p_lik <- ggplot(results, aes(x = rho, y = logLik)) +

  # Fitted profile likelihood curve
  geom_line(linewidth = 0.9, color = "black") +
  geom_point(size = 2.5, color = "black") +

  # Vertical line at the protocol assumption (rho = 0.5)
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    color      = "grey50",
    linewidth  = 0.7
  ) +
  annotate("text", x = 0.5, y = min(results$logLik, na.rm = TRUE),
           label = "rho = 0.5\n(protocol)", hjust = -0.1, vjust = 0,
           size = 3, color = "grey40") +

  # Vertical line at the best rho (if different from 0.5)
  {if (best_rho != 0.5)
    geom_vline(
      xintercept = best_rho,
      linetype   = "solid",
      color      = "#d6604d",
      linewidth  = 0.7
    )
  } +
  {if (best_rho != 0.5)
    annotate("text", x = best_rho, y = min(results$logLik, na.rm = TRUE),
             label = paste0("rho = ", best_rho, "\n(best fit)"),
             hjust = 1.1, vjust = 0, size = 3, color = "#d6604d")
  } +

  labs(
    x        = "Within-study correlation (rho)",
    y        = "REML log-likelihood",
    title    = "Profile REML likelihood over rho",
    subtitle = paste0(
      "CHE model (intercept only). Best fit at rho = ", best_rho,
      "; protocol uses rho = 0.5."
    )
  ) +
  scale_x_continuous(breaks = seq(0, 0.95, by = 0.1)) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(size = 9))

ggsave(
  filename = "Figure_RhoLikelihood.png",
  plot     = p_lik,
  width    = 7,
  height   = 5,
  dpi      = 300
)
cat("\nFigure_RhoLikelihood.png saved.\n")

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

  # Vertical line at rho = 0.5
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    color      = "grey50",
    linewidth  = 0.7
  ) +
  annotate("text", x = 0.5,
           y = min(results_plot$ci_lo, na.rm = TRUE),
           label = "rho = 0.5", hjust = -0.1, vjust = 0,
           size = 3, color = "grey40") +

  # Vertical line at best rho (if different from 0.5)
  {if (best_rho != 0.5)
    geom_vline(
      xintercept = best_rho,
      linetype   = "solid",
      color      = "#d6604d",
      linewidth  = 0.7
    )
  } +
  {if (best_rho != 0.5)
    annotate("text", x = best_rho,
             y = min(results_plot$ci_lo, na.rm = TRUE),
             label = paste0("rho = ", best_rho),
             hjust = 1.1, vjust = 0, size = 3, color = "#d6604d")
  } +

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
  filename = "Figure_RhoIntercept.png",
  plot     = p_int,
  width    = 7,
  height   = 5,
  dpi      = 300
)
cat("Figure_RhoIntercept.png saved.\n")

# =============================================================================
# Results table
# =============================================================================

tbl <- results %>%
  mutate(
    best_fit = rho == best_rho,
    across(c(logLik, intercept, se, tau, omega), ~ round(.x, 4))
  ) %>%
  rename(
    `Rho`           = rho,
    `REML logLik`   = logLik,
    `Intercept`     = intercept,
    `SE (CR2)`      = se,
    `Tau`           = tau,
    `Omega`         = omega,
    `Best fit`      = best_fit
  )

write_xlsx(
  list("Rho sensitivity" = tbl),
  path = "Table_RhoSensitivity.xlsx"
)
cat("Table_RhoSensitivity.xlsx saved.\n")

cat("\nDone. Summary:\n")
cat(sprintf("  Best rho:             %.2f\n", best_rho))
cat(sprintf("  Intercept at rho=0.5: %.4f (SE = %.4f)\n",
            results$intercept[results$rho == 0.50],
            results$se[results$rho == 0.50]))
cat(sprintf("  Intercept at best rho:%.4f (SE = %.4f)\n",
            results$intercept[best_idx],
            results$se[best_idx]))

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
              