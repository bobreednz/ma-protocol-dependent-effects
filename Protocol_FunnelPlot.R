# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("here", "metafor", "haven")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)      # resolves file paths relative to the project root (.Rproj)
library(metafor)  # rma() and funnel()

# ============================================================
# Load data
# ============================================================

DT <- readRDS(here("SCData_processed.rds"))

# ============================================================
# Fit a simple random effects model
#
# We use rma() here (not rma.mv()) because funnel plots are a
# descriptive/diagnostic tool that does not require the full
# CHE specification. rma() fits a standard RE model with one
# variance component, treating each effect size as independent.
# In a dataset with multiple effect sizes per study this
# slightly underestimates uncertainty, but it is conventional
# for funnel plot construction.
#
# The RE estimate serves as the center of the funnel.
# ============================================================

model_RE <- rma(yi = z, sei = sez, data = DT)

cat(sprintf("RE estimate used as funnel center: %.3f (SE = %.3f)\n",
            model_RE$beta, model_RE$se))

# ============================================================
# FIGURE 2: Standard funnel plot
#
# Standard error (sez) is on the y-axis with 0 at the top and
# increasing downward -- more precise estimates (smaller SE)
# sit near the top, less precise ones near the bottom. The
# shaded region is the 95% pseudo-confidence interval around
# the estimated mean. Under no publication bias and no
# heterogeneity, points should be roughly symmetrically
# distributed within the funnel. The rightward asymmetry
# visible here -- a concentration of points to the right of
# the funnel center, especially among less precise studies --
# is a visual signal of publication bias, motivating the
# formal PET-PEESE tests in Section 5.2.
# ============================================================

png(here("Figure_FunnelPlot.png"),
    width = 6, height = 5, units = "in", res = 300)

par(mar = c(4, 4, 2, 1.5))

# cex = 0.75 draws the points 25% smaller than the default, so that the
# dense cluster of precise estimates near the top of the funnel is easier
# to see.
metafor::funnel(model_RE,
       cex  = 0.75,
       xlab = "Fisher's z",
       ylab = "Standard Error")

dev.off()

cat("Figure 2 saved to Figure_FunnelPlot.png\n")

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
