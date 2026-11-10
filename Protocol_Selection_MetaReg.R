# =============================================================================
# Protocol_Selection_MetaReg.R
# =============================================================================
# Purpose: the step-function selection model with the 18 moderators, and its
# best-practice predictions (Section 7.2).
#
# Which selection model? Section 6.2 decides between the 3PSM and the 4PSM.
# Set psm_model below to the one it favors:
#   "3PSM": one threshold at one-sided p = 0.025
#   "4PSM": thresholds at one-sided p = 0.025 and 0.5 (used in this example,
#           because the diagnostics show selection on sign)
#
# The moderators enter the mean of the distribution of effects BEFORE
# selection. The selection weights (lambdas) and the heterogeneity (tau) are
# common to all estimates, as in Section 5.3.
#
# Three fits, identical except for where the moderators are centered:
#   "mean": moderators centered at their sample means. The intercept is the
#           corrected mean effect for a study with average characteristics
#           (comparable with the CHE-PET-PEESE intercept in Section 7.1).
#           The lambdas and the moderator coefficients come from this fit.
#   "bp1":  moderators centered at the best-practice values, non-OECD/Europe.
#   "bp2":  the same, OECD/Europe.
#           In the two best-practice fits the intercept IS the best-practice
#           prediction, so its bootstrap interval comes directly from the fit.
#           (Centering changes only the intercept, so the three fits share the
#           same moderator coefficients and selection weights.)
#
# Continuous moderators (PubYear, NumberSCVars) are also divided by their
# standard deviation before estimation. This only helps the optimizer: the
# table converts their coefficients back to effects per year and per
# additional social-capital variable.
#
# Intervals: percentile intervals from a two-stage cluster bootstrap
# (R = 1999), as in Section 5.3; no large-sample intervals are reported.
# A bootstrap replication fails when the resampled data no longer identify
# every coefficient, for example when no study with a rare moderator is
# drawn. The intervals use the replications that converged; the table
# reports how many there were.
#
# Prediction interval: the range of TRUE effects among best-practice
# estimates. For each bootstrap replication, true effects are drawn from a
# normal distribution with that replication's best-practice mean and tau^2;
# the 2.5th and 97.5th percentiles of all draws give the 95% interval. It
# does not include the sampling error of a new estimate.
#
# Estimand: like the RE model, the selection model treats each estimate as
# marginally independent, so its mean is an average across estimates rather
# than the CHE mean (see the note in Section 5.3).
#
# Output: Table_Selection_MetaReg.xlsx
#   "Summary"               corrected mean at average characteristics,
#                           selection weights, tau, replications used
#   "Moderators"            moderator coefficients with bootstrap intervals,
#                           and the number of studies with each binary
#                           characteristic (a rare one makes replications fail)
#   "BP1 - Non-OECD Europe" best-practice prediction, Fisher's z and PCC
#   "BP2 - OECD Europe"     the same for OECD/Europe
#   "Notes"                 the best-practice values, one note per scenario
#   "Values"                unformatted numbers, for use in Section 8
# Fitted models: sel_4psm_reg18_mean.rds, sel_4psm_reg18_bp1.rds,
#   sel_4psm_reg18_bp2.rds (sel_3psm_... for the 3PSM).
#
# Run time: with refit <- FALSE (the default), a few seconds; the script
# loads the saved fits. With refit <- TRUE, about 15 to 20 minutes with four
# cores (three bootstrap fits of about 5 minutes each).
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

# ---- Settings ------------------------------------------------------------------
psm_model <- "4PSM"   # "4PSM" or "3PSM": the model favored in Section 6.2
refit     <- FALSE    # TRUE refits the three models (slow); FALSE loads saved fits

steps <- switch(psm_model,
                "3PSM" = 0.025,
                "4PSM" = c(0.025, 0.5),
                stop('psm_model must be "3PSM" or "4PSM".'))

# ---- Packages ----------------------------------------------------------------
# metaselection is on CRAN (version 0.3.0 or later). The results reported in
# the protocol were produced with development version 0.1.5.
packages <- c("here", "metaselection", "future", "progressr", "writexl")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)      # resolves file paths relative to the project root (.Rproj)
library(writexl)   # write_xlsx()
# metaselection, future and progressr are loaded below, only if refitting

# ---- Data --------------------------------------------------------------------
dat <- readRDS(here("SCData_processed.rds"))
n_obs     <- nrow(dat)
n_studies <- length(unique(dat$newid))
cat(sprintf("Dataset: %d observations from %d studies.\n", n_obs, n_studies))

# ---- Moderators ---------------------------------------------------------------
# Same 18 moderators, in the same order, as in Section 7.1
moderators <- c(
  "SC1_Cognitive", "SC1_Structural", "DV_GrowthRate", "PubYear", "Published",
  "LaggedDV", "LaggedSC", "NumberSCVars", "Endog_IV", "Endog_FE",
  "CityLevel", "RegionLevel", "CountryLevel", "PanelData",
  "Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia"
)

# Labels used in the table (same wording as Table 1 and Section 7.1)
mod_labels <- c(
  SC1_Cognitive  = "SC: Cognitive",
  SC1_Structural = "SC: Structural",
  DV_GrowthRate  = "Outcome: Growth rate",
  PubYear        = "Publication year",
  Published      = "Published",
  LaggedDV       = "Lagged dependent variable",
  LaggedSC       = "Lagged social capital",
  NumberSCVars   = "Number of SC variables",
  Endog_IV       = "Endogeneity: IV",
  Endog_FE       = "Endogeneity: FE",
  CityLevel      = "City-level",
  RegionLevel    = "Region-level",
  CountryLevel   = "Country-level",
  PanelData      = "Panel data",
  Reg_OECDEurope = "Region: OECD/Europe",
  Reg_US         = "Region: US",
  Reg_Africa     = "Region: Africa",
  Reg_Asia       = "Region: Asia"
)

# Continuous moderators are divided by their standard deviation
continuous <- c("PubYear", "NumberSCVars")
sds <- sapply(continuous, function(v) sd(as.numeric(dat[[v]])))

# Best-practice values (identical to Protocol_CHE_PETPEESE_BestPractice.R). The two scenarios
# differ only in region: BP1 = non-OECD/Europe, BP2 = OECD/Europe.
bp_values <- function(oecd) c(
  SC1_Cognitive = 1, SC1_Structural = 0, DV_GrowthRate = 1,
  PubYear = mean(as.numeric(dat$PubYear)), Published = 1,
  LaggedDV = 0, LaggedSC = 0, NumberSCVars = 1, Endog_IV = 1, Endog_FE = 1,
  CityLevel = 0, RegionLevel = 0, CountryLevel = 1, PanelData = 1,
  Reg_OECDEurope = oecd, Reg_US = 0, Reg_Africa = 0, Reg_Asia = 0)
mean_values <- sapply(moderators, function(v) mean(as.numeric(dat[[v]])))

# Center each moderator at the chosen values: x_c = (x - centre) / scale,
# where scale is the SD for continuous moderators and 1 for binary ones. The
# intercept is then the corrected effect at the centering values.
make_data <- function(centre) {
  d <- dat
  for (v in moderators) {
    scl <- if (v %in% continuous) sds[[v]] else 1
    d[[paste0(v, "_c")]] <- (as.numeric(d[[v]]) - centre[[v]]) / scl
  }
  d
}
centres <- list(mean = mean_values, bp1 = bp_values(0), bp2 = bp_values(1))

fit_file <- function(which) here(sprintf("sel_%s_reg18_%s.rds", tolower(psm_model), which))

# ---- Check for rare moderators before estimating -------------------------------
# A bootstrap replication fails when a moderator takes the same value in every
# resampled estimate. This happens when none of the studies with (or without)
# a characteristic is drawn. If m of J studies have the characteristic, the
# chance that a resample of J studies contains none of them is (1 - m/J)^J,
# roughly exp(-m): about 37% for m = 1, 13% for m = 2 and 5% for m = 3.
# The table printed below reports this for each binary moderator, so the
# researcher can see in advance which moderators will cause replications to
# fail. The number of studies also appears in the moderator table.
binary <- setdiff(moderators, continuous)
study_counts <- do.call(rbind, lapply(binary, function(v) {
  x <- as.numeric(dat[[v]])
  with_1 <- sum(tapply(x, dat$newid, max) == 1)   # studies with at least one estimate = 1
  with_0 <- sum(tapply(x, dat$newid, min) == 0)   # studies with at least one estimate = 0
  # Share of resamples in which the moderator is constant (approximate; ignores
  # the second, within-study stage of the bootstrap)
  p_const <- (1 - with_1 / n_studies)^n_studies + (1 - with_0 / n_studies)^n_studies
  data.frame(Moderator = unname(mod_labels[v]),
             `Studies with characteristic` = with_1,
             `Studies without` = with_0,
             `Resamples where constant (%)` = round(100 * p_const, 1),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
cat("\nStudy counts for the binary moderators:\n")
print(study_counts, row.names = FALSE)

# ---- Fit or load the three models ----------------------------------------------
if (refit) {
  library(metaselection)  # selection_model() and define_priors()
  library(future)         # parallel processing for the bootstrap
  library(progressr)      # progress bar during the bootstrap

  # Same weak priors as Section 5.3 (Pustejovsky, Citkowicz, and Joshi, 2026)
  prior_spec <- define_priors(
    beta_mean = 0, beta_precision = 1 / 16, beta_L = 4,
    tau_mode = 0.2, tau_alpha = 1,
    lambda_mode = 0.5, lambda_precision = 1 / 54, lambda_L = 4
  )
  mod_formula <- as.formula(paste("~", paste(paste0(moderators, "_c"), collapse = " + ")))

  handlers(global = TRUE)
  plan(multisession, workers = 4)   # adjust to the number of cores available

  for (which in names(centres)) {
    cat(sprintf("\nFitting the %s meta-regression, moderators centered at %s...\n",
                psm_model, which))
    set.seed(12345)   # the same seed for each fit, as in Section 5.3
    fit <- selection_model(
      data = make_data(centres[[which]]), yi = z, sei = sez, cluster = newid,
      selection_type = "step",
      steps          = steps,
      mean_mods      = mod_formula,                 # moderators enter the mean only
      estimator      = "CML",                       # with the priors, this is the PML estimator
      priors         = prior_spec,
      bootstrap      = "two-stage",                 # resamples studies, then estimates
      CI_type        = c("large-sample", "percentile"),
      R              = 1999
    )
    saveRDS(fit, fit_file(which))
  }
  plan(sequential)   # close the parallel workers
}

fits <- lapply(names(centres), function(which) {
  f <- fit_file(which)
  if (!file.exists(f)) stop("Saved fit not found: ", f, ". Set refit <- TRUE.")
  readRDS(f)
})
names(fits) <- names(centres)

# Guard against a mismatch between psm_model and the saved fits
if (!isTRUE(all.equal(as.numeric(fits$mean$steps), steps)))
  stop("The saved fits do not use the thresholds of ", psm_model, ".")

# ---- Extract estimates -----------------------------------------------------------
# Each fit stores its estimates in fit$est, one row per parameter: beta_* for
# the mean coefficients, gamma = log(tau^2), and zeta1, zeta2 = log(lambda1),
# log(lambda2). Columns Est, percentile_lower and percentile_upper hold the
# estimate and the bootstrap percentile interval. Because the transformations
# are monotonic, exponentiating the interval ends gives the intervals for
# tau^2 and the lambdas.
get_row <- function(fit, name) {
  r <- fit$est[fit$est$param == name, ]
  if (nrow(r) == 0) NULL else r
}
estimate <- function(fit, name, transform = identity) {
  r <- get_row(fit, name)
  if (is.null(r)) return(c(est = NA, lo = NA, hi = NA))
  c(est = transform(r$Est), lo = transform(r$percentile_lower),
    hi = transform(r$percentile_upper))
}

# Number of bootstrap replications that converged (out of the 1,999 requested)
n_boot <- function(fit) {
  if ("bootstraps" %in% names(fit$est)) fit$est$bootstraps[1] else fit$R
}

# Bootstrap prediction interval for the true effects at the centering values
# (see the header). Pairs each replication's intercept with its tau^2.
predict_interval <- function(fit, n_draw = 200, seed = 20261016) {
  reps <- fit$bootstrap_reps
  b0   <- reps$Est[reps$param == "beta_(Intercept)"]
  g    <- reps$Est[reps$param == "gamma"]
  set.seed(seed)
  draws <- rnorm(length(b0) * n_draw, mean = rep(b0, each = n_draw),
                 sd = rep(sqrt(exp(g)), each = n_draw))
  unname(quantile(draws, c(0.025, 0.975)))
}

# ---- Formatting helpers ---------------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "--", formatC(x, digits = 3, format = "f"))
f4 <- function(x) ifelse(is.na(x), "--", formatC(x, digits = 4, format = "f"))
ci <- function(lo, hi, f = f3) ifelse(is.na(lo), "--", paste0("[", f(lo), ", ", f(hi), "]"))

# ---- Sheet 1: Summary ------------------------------------------------------------
fm  <- fits$mean
mu  <- estimate(fm, "beta_(Intercept)")
l1  <- estimate(fm, "zeta1", exp)
l2  <- estimate(fm, "zeta2", exp)                    # NA for the 3PSM
tau <- estimate(fm, "gamma", function(x) sqrt(exp(x)))

summary_tbl <- data.frame(
  Parameter = c("Corrected mean", "95% CI (bootstrap)",
                "lambda1", "95% CI (bootstrap)",
                "lambda2", "95% CI (bootstrap)",
                "tau (total)", "Bootstrap replications used",
                "Observations", "Studies"),
  Value = c(f3(mu["est"]), ci(mu["lo"], mu["hi"]),
            f3(l1["est"]), ci(l1["lo"], l1["hi"]),
            f3(l2["est"]), ci(l2["lo"], l2["hi"]),
            f3(tau["est"]),
            sprintf("%s of 1,999", format(n_boot(fm), big.mark = ",")),
            n_obs, n_studies),
  stringsAsFactors = FALSE
)
names(summary_tbl)[2] <- psm_model

# ---- Sheet 2: Moderator coefficients ----------------------------------------------
# Continuous moderators: divide by the SD to get effects per year (PubYear)
# and per additional social-capital variable (NumberSCVars). Dividing the
# interval ends by the same positive number gives the rescaled interval.
mod_est <- t(sapply(moderators, function(v) {
  e <- estimate(fm, paste0("beta_", v, "_c"))
  if (v %in% continuous) e / sds[[v]] else e
}))
small_digits <- moderators %in% continuous   # per-unit effects are small: show 4 decimals
# Number of studies with each binary characteristic (blank for continuous
# moderators), taken from the rare-moderator check above
n_with <- setNames(study_counts$`Studies with characteristic`, binary)
moderator_tbl <- data.frame(
  Moderator            = unname(mod_labels[moderators]),
  Studies              = ifelse(moderators %in% binary, as.character(n_with[moderators]), "--"),
  Coefficient          = ifelse(small_digits, f4(mod_est[, "est"]), f3(mod_est[, "est"])),
  `95% CI (bootstrap)` = ifelse(small_digits, ci(mod_est[, "lo"], mod_est[, "hi"], f4),
                                ci(mod_est[, "lo"], mod_est[, "hi"])),
  check.names = FALSE, stringsAsFactors = FALSE
)

# ---- Sheets 3 and 4: best-practice predictions ------------------------------------
bp_panel <- function(fit) {
  e  <- estimate(fit, "beta_(Intercept)")
  pi <- predict_interval(fit)
  list(
    display = data.frame(
      `Effect Size`     = c("Fisher's z", "PCC"),
      `Mean Prediction` = c(f3(e["est"]), f3(tanh(e["est"]))),
      `95% CI`          = c(ci(e["lo"], e["hi"]), ci(tanh(e["lo"]), tanh(e["hi"]))),
      `95% PI`          = c(ci(pi[1], pi[2]), ci(tanh(pi[1]), tanh(pi[2]))),
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    values = c(e, pi_lo = pi[1], pi_hi = pi[2], n_boot = n_boot(fit))
  )
}
bp1 <- bp_panel(fits$bp1)
bp2 <- bp_panel(fits$bp2)

# Notes listing every best-practice value, as in Protocol_CHE_PETPEESE_BestPractice.R
bp_note <- function(oecd, fit) {
  paste0(
    "\"Best Practice\" predictions are derived from the ", psm_model,
    " meta-regression with 18 moderators (Section 7.2). ",
    "Predictions are evaluated at SC1_Cognitive = 1, SC1_Structural = 0, ",
    "SC1_Other = 0 (reference), DV_GrowthRate = 1, DV_GDPLevel = 0 (reference), PubYear = ",
    round(mean(as.numeric(dat$PubYear)), 1), " (sample mean), Published = 1, LaggedDV = 0, ",
    "LaggedSC = 0, NumberSCVars = 1, Endog_IV = 1, Endog_FE = 1, CityLevel = 0, RegionLevel = 0, ",
    "CountryLevel = 1, PanelData = 1, Reg_OECDEurope = ", oecd,
    ", Reg_US = 0, Reg_Africa = 0, Reg_Asia = 0. ",
    "The prediction is the mean of the distribution of effects before selection. ",
    "Intervals are percentile intervals from a two-stage cluster bootstrap, based on the ",
    format(n_boot(fit), big.mark = ","), " of 1,999 replications that converged. ",
    "The prediction interval is the range of true effects (bootstrap draws from a normal ",
    "distribution with the predicted mean and tau^2); it does not include sampling error."
  )
}
notes_sheet <- data.frame(
  Scenario = c("BP1 - Non-OECD/Europe", "BP2 - OECD/Europe"),
  Note     = c(bp_note(0, fits$bp1), bp_note(1, fits$bp2)),
  stringsAsFactors = FALSE
)

# ---- Sheet 6: unformatted values ----------------------------------------------------
value_row <- function(quantity, e, pi = c(NA, NA)) {
  data.frame(Model = psm_model, Quantity = quantity,
             Estimate = unname(e["est"]), Lower = unname(e["lo"]), Upper = unname(e["hi"]),
             PI_lower = pi[1], PI_upper = pi[2], stringsAsFactors = FALSE)
}
values <- rbind(
  value_row("Corrected mean", mu),
  value_row("lambda1", l1),
  value_row("lambda2", l2),
  value_row("tau", tau),
  do.call(rbind, lapply(seq_along(moderators), function(i)
    value_row(unname(mod_labels[moderators[i]]), mod_est[i, ]))),
  value_row("BP1 (non-OECD/Europe)", bp1$values[c("est", "lo", "hi")],
            unname(bp1$values[c("pi_lo", "pi_hi")])),
  value_row("BP2 (OECD/Europe)", bp2$values[c("est", "lo", "hi")],
            unname(bp2$values[c("pi_lo", "pi_hi")]))
)
values$Replications <- c(rep(n_boot(fm), 4 + length(moderators)),
                         bp1$values["n_boot"], bp2$values["n_boot"])

# ---- Save -------------------------------------------------------------------------
write_xlsx(list(Summary                 = summary_tbl,
                Moderators              = moderator_tbl,
                "BP1 - Non-OECD Europe" = bp1$display,
                "BP2 - OECD Europe"     = bp2$display,
                Notes                   = notes_sheet,
                Values                  = values),
           path = here("Table_Selection_MetaReg.xlsx"))

print(summary_tbl, row.names = FALSE)
print(moderator_tbl, row.names = FALSE)
cat("\nBP1 (non-OECD/Europe):\n"); print(bp1$display, row.names = FALSE)
cat("\nBP2 (OECD/Europe):\n");     print(bp2$display, row.names = FALSE)
cat(sprintf("\nSaved Table_Selection_MetaReg.xlsx. Run time: %.1f seconds.\n",
            (proc.time() - start_time)["elapsed"]))
