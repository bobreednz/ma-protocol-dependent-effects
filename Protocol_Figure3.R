library(here)   # resolves file paths relative to the project root (.Rproj)
library(haven)
library(RoBMA)

# ---- LOAD DATA ----
dat <- readRDS(here("SCData_processed.rds"))

# ---- TRANSFORM TO FISHER'S Z ----
# atanh() applies the Fisher's z transformation to the partial correlations.
# Variance of Fisher's z for a partial correlation is 1 / (df - 1).
dat$yi <- atanh(dat$pcc)
dat$vi <- 1 / (dat$df - 1)

# ---- FIT INTERCEPT-ONLY MULTILEVEL RoBMA ----
# cluster = newid specifies effect sizes nested within studies.
# measure = "ZCOR" tells RoBMA the input is Fisher's z transformed correlations.
#
# ---- ADAPTING THIS SCRIPT FOR A DIFFERENT EFFECT SIZE ----
# Fisher's z is one of six standardized effect-size types RoBMA supports
# directly via the measure argument, with no further setup: "SMD", "ZCOR",
# "RR", "OR", "HR", and "IRR". Each has a known unit information standard
# deviation (UISD) built into the package, which is what lets measure =
# "ZCOR" work here without any additional scaling step.
#
# A nonstandard effect size without a known UISD -- a raw (unstandardized)
# regression coefficient being the most common example a researcher adapting
# this protocol is likely to encounter -- needs measure = "GEN" instead, plus
# one of the following two approaches (per direct correspondence with
# Frantisek Bartos, RoBMA's developer, 2026-07-28, and confirmed against
# "RoBMA documentation.pdf" in this folder):
#
#   1. If sample sizes are available for every effect size, pass them as an
#      additional ni argument alongside the usual yi/vi:
#
#        fit <- RoBMA(
#          yi      = yi,       # raw regression coefficients
#          vi      = vi,       # their squared standard errors
#          ni      = ni,       # sample size for each estimated coefficient
#          measure = "GEN",
#          mods    = ~ ...,
#          cluster = newid,
#          ...
#        )
#
#      RoBMA uses ni (together with vi) to estimate the UISD automatically
#      and scale its default priors accordingly -- no other change is needed.
#
#   2. If sample sizes are only available for some of the effect sizes,
#      compute the UISD once from that subset via
#      estimate_unit_information_sd() and pass the resulting number in
#      through prior_unit_information_sd, instead of supplying ni for every
#      observation:
#
#        uisd <- estimate_unit_information_sd(sei = sei_subset, ni = ni_subset)
#        fit  <- RoBMA(
#          yi = yi, vi = vi, measure = "GEN",
#          prior_unit_information_sd = uisd,
#          mods = ~ ..., cluster = newid, ...
#        )
#
#      estimate_unit_information_sd() requires complete (non-missing) sei/ni
#      vectors, so it is run only on the subset of studies where both are
#      actually available; the single UISD value it returns is then reused
#      for the full model.
#
# Two further points:
#
#   - This only affects RoBMA (this script, Protocol_Table10.R, and
#     Protocol_Table10_FYI.R). The CHE/FAT-PET machinery in Tables 4-7 is
#     built on metafor, which only needs yi and vi and is indifferent to
#     what scale the effect size is on.
#   - The atanh()/tanh() transform-and-back-transform steps used throughout
#     this protocol (see TRANSFORM TO FISHER'S Z above, and the PCC
#     back-transform in the Table 7/11 scripts) are specific to correlational
#     effect sizes and should simply be omitted for a raw regression
#     coefficient or any other non-correlational measure.

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
  saveRDS(fit, here("RoBMA_intercept_only.rds"))
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
