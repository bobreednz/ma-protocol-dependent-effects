library(haven)
library(RoBMA)

setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# ---- LOAD DATA ----
dat <- readRDS("SCData_processed.rds")

# ---- TRANSFORM TO FISHER'S Z ----
# atanh() applies the Fisher's z transformation to the partial correlations.
# Variance of Fisher's z for a partial correlation is 1 / (df - 1).
dat$yi <- atanh(dat$pcc)
dat$vi <- 1 / (dat$df - 1)

# ---- FIT INTERCEPT-ONLY MULTILEVEL RoBMA ----
# cluster = newid specifies effect sizes nested within studies.
# measure = "ZCOR" tells RoBMA the input is Fisher's z transformed correlations.
cat("Fitting intercept-only multilevel RoBMA...\n")
time_fit <- system.time({
  fit <- RoBMA(
    yi = yi, vi = vi, measure = "ZCOR",
    cluster = newid,
    sample = 20000, burnin = 10000, adapt = 10000,
    thin = 5, parallel = TRUE, seed = 1,
    data = dat
  )
})
cat("Model fitting time (seconds):\n")
print(time_fit)

# ---- SAVE MODEL OBJECT ----
# Save so Protocol_Table8.R can load without refitting.
cat("Saving model object...\n")
time_save <- system.time({
  saveRDS(fit, "RoBMA_intercept_only.rds")
})
cat("Save time (seconds):\n")
print(time_save)

# ---- SUMMARY ----
cat("Model summary:\n")
time_summary <- system.time({
  print(summary(fit, include_mcmc_diagnostics = FALSE))
})
cat("Summary time (seconds):\n")
print(time_summary)

# ---- ZPLOT: INTERACTIVE ----
cat("Producing zplot...\n")
time_zplot <- system.time({
  par(mar = c(4, 4, 0, 0))
  zplot(fit, by.hist = 0.25, plot_extrapolation = FALSE, from = -4, to = 8)
})
cat("Zplot time (seconds):\n")
print(time_zplot)

# ---- ZPLOT: SAVE AS PNG ----
# dev.copy() replays the plot already drawn on the device above, avoiding a
# second (expensive, ~10 minute) zplot() computation.
cat("Saving zplot as PNG...\n")
time_png <- system.time({
  dev.copy(png, "Figure3_Zplot.png", width = 800, height = 600, res = 120)
  dev.off()
})
cat("PNG save time (seconds):\n")
print(time_png)

cat("Protocol_Figure3.R complete.\n")
