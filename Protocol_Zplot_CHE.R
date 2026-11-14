# =============================================================================
# Protocol_Zplot_CHE.R
# =============================================================================
# Purpose: z-plot for the intercept-only CHE RoBMA model (Section 6).
#
# The z-plot (Bartos and Schimmack, 2026) shows the histogram of the observed
# z-statistics (each estimate divided by its standard error) together with the
# distribution of z-statistics implied by the fitted RoBMA model. It serves
# two purposes in the protocol:
#   1. Diagnosing selection. A sharp change in the histogram at z = 0 points
#      to selection on sign; a sharp change at the significance thresholds
#      (z = 1.645 and 1.96) points to selection on significance.
#   2. Checking fit. Where the histogram departs from the fitted curve, the
#      model does not reproduce that part of the observed distribution.
#
# The script loads the fit saved by Protocol_RoBMA_CHE_Intercept.R; it does
# not refit the model. zplot() takes about 5 minutes here.
#
# Output: Figure_Zplot_CHE.png
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

library(here)    # resolves file paths relative to the project root (.Rproj)
library(RoBMA)   # zplot(); requires RoBMA 4.1.1 or later for RoBMA.mv fits

# The CHE models need RoBMA 4.1.1 or later, the first version with RoBMA.mv().
if (packageVersion("RoBMA") < "4.1.1") {
  stop("This script needs RoBMA 4.1.1 or later. ",
       "Run Protocol_Install.R, or see Section 5.4 of the protocol.")
}

fit_rob <- readRDS(here("fit_RoBMA_CHE_intercept.rds"))

# Draw the z-plot on screen. by.hist sets the histogram bin width; from/to set
# the range of z shown (a few estimates lie outside it); the extrapolation of
# the fitted curve beyond the observed data is switched off.
par(mar = c(4, 4, 0, 0))
zplot(fit_rob, by.hist = 0.25, plot_extrapolation = FALSE, from = -4, to = 8)

# Copy the plot already on screen to a PNG. This avoids recomputing it.
dev.copy(png, here("Figure_Zplot_CHE.png"), width = 800, height = 600, res = 120)
dev.off()
cat("Saved Figure_Zplot_CHE.png\n")

elapsed <- proc.time() - start_time
cat(sprintf("Total run time: %.1f minutes.\n", elapsed["elapsed"] / 60))
