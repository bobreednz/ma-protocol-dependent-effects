# =============================================================================
# Protocol_RoBMA_CHE_Intercept.R
# =============================================================================
# Purpose: intercept-only CHE RoBMA (Section 5.4), with an uncorrected Bayesian
# benchmark.
#
#   RoBMA.mv: Robust Bayesian model-averaged meta-analysis with the CHE
#             dependence structure. It averages over models with and without
#             an effect, with and without heterogeneity, and with and without
#             publication bias. The publication-bias models are six
#             step-function selection models plus PET and PEESE.
#   BMA.mv:   the same model with the publication-bias component switched off.
#             It is the Bayesian counterpart of the uncorrected CHE mean in the
#             five-estimators table. Comparing BMA.mv with RoBMA.mv shows how
#             much of the change in the mean is due to the bias correction.
#
# Both use the CHE structure of the frequentist models: correlated sampling
# errors (V, built by vcalc() with rho = 0.5) plus study-level and
# estimate-level random effects (random = ~ 1 | newid/obs).
#
# SOFTWARE: RoBMA.mv() and BMA.mv() require RoBMA 4.1, which at the time of
# writing is a development version (the CRAN release is 4.0.0). Install it,
# then restart R:
#   remotes::install_github("FBartos/BayesTools@0.3.1-add-random")
#   remotes::install_github("FBartos/RoBMA@dev-4-1-0")
# RoBMA also needs JAGS (see Software Setup in the Introduction).
#
# SAVED FITS: fitting takes about 4.5 hours for RoBMA.mv and a few minutes for
# BMA.mv. By default the script loads the saved fits. Set refit <- TRUE to
# re-estimate them (the fits are then saved again).
#
# Outputs:
#   fit_RoBMA_CHE_intercept.rds, fit_BMA_CHE_intercept.rds   fitted models
#   Table_RoBMA_CHE_Intercept.xlsx
#     "Components":       prior and posterior probabilities, inclusion Bayes factors
#     "Estimates":        mean effect and heterogeneity, RoBMA and BMA
#     "SelectionWeights": RoBMA publication-bias weights
#     "Values":           numeric values used by Protocol_Diagnostics.R
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

refit <- FALSE              # TRUE re-estimates both models (about 4.5 hours)

# ---- Packages ----------------------------------------------------------------
library(here)      # resolves file paths relative to the project root (.Rproj)
library(metafor)   # vcalc() for the sampling covariance matrix
library(RoBMA)     # RoBMA.mv() and BMA.mv(); must be version 4.1 or later
library(writexl)   # write_xlsx()

# The CHE models need RoBMA 4.1 or later: RoBMA.mv() is not in the CRAN version (4.0.0).
if (packageVersion("RoBMA") < "4.1") {
  stop("This script needs the development version of RoBMA (4.1 or later). ",
       "Run Protocol_Install.R, or see Section 5.4 of the protocol.")
}

cat("RoBMA version:", as.character(packageVersion("RoBMA")), "\n")

rob_file <- here("fit_RoBMA_CHE_intercept.rds")
bma_file <- here("fit_BMA_CHE_intercept.rds")

# ---- Fit (only if requested or if the saved fits are missing) -------------
if (refit || !file.exists(rob_file) || !file.exists(bma_file)) {

  dat <- readRDS(here("SCData_processed.rds"))

  # Fisher's z. For other effect sizes, see the note at the end of this script.
  dat$yi <- atanh(dat$pcc)

  # Correlated sampling errors (the "C" in CHE), as in the frequentist models
  V_mat <- vcalc(vi = sez^2, cluster = newid, obs = obs, rho = 0.5, data = dat)

  # MCMC settings. thin = 2 runs the chains twice as long and keeps every
  # second draw, which raises the effective sample size of the mean.
  mcmc <- list(sample = 10000, burnin = 10000, adapt = 10000, thin = 2)

  cat("Fitting intercept-only CHE RoBMA (about 4.5 hours)...\n")
  t1 <- system.time({
    fit_rob <- RoBMA.mv(
      yi = yi, V = V_mat, measure = "ZCOR",   # input is Fisher's z
      random = ~ 1 | newid/obs,               # study and estimate random effects
      sample = mcmc$sample, burnin = mcmc$burnin, adapt = mcmc$adapt,
      thin = mcmc$thin, parallel = TRUE, seed = 1, data = dat)
  })
  cat("RoBMA fitting time (hours):", round(t1["elapsed"] / 3600, 2), "\n")
  saveRDS(fit_rob, rob_file)

  cat("Fitting intercept-only CHE BMA, no bias component (a few minutes)...\n")
  t2 <- system.time({
    fit_bma <- BMA.mv(
      yi = yi, V = V_mat, measure = "ZCOR",
      random = ~ 1 | newid/obs,
      sample = mcmc$sample, burnin = mcmc$burnin, adapt = mcmc$adapt,
      thin = mcmc$thin, parallel = TRUE, seed = 1, data = dat)
  })
  cat("BMA fitting time (minutes):", round(t2["elapsed"] / 60, 1), "\n")
  saveRDS(fit_bma, bma_file)

} else {
  cat("Loading saved fits (set refit <- TRUE to re-estimate)...\n")
  fit_rob <- readRDS(rob_file)
  fit_bma <- readRDS(bma_file)
}

# ---- Summaries ---------------------------------------------------------------
# conditional = TRUE adds estimates conditional on the presence of an effect
# to the usual model-averaged estimates. Printing the summaries also shows the
# MCMC diagnostics (R-hat and effective sample size) for every parameter.
s_rob <- summary(fit_rob, conditional = TRUE)
s_bma <- summary(fit_bma, conditional = TRUE)
print(s_rob)
print(s_bma)

# The summary object is a list of tables. Its element names are printed here
# so that, if a future version of RoBMA renames them, the fix is easy to find.
cat("\nElements of the RoBMA summary object:", paste(names(s_rob), collapse = ", "), "\n")

# ---- Helpers for reading the summary tables ----------------------------------
# In RoBMA 4.1, a RoBMA.mv summary stores (element names checked 2026-10-14):
#   inclusion_components  Effect and Publication Bias inclusion
#   inclusion_random      heterogeneity (tau_total) inclusion
#   estimates_mods        the mean effect, mu (model-averaged)
#   estimates_mods_conditional   mu conditional on an effect
#   estimates_random      tau_total and the split of heterogeneity
#   estimates_bias        selection weights (omega) and PET/PEESE coefficients
# get_tab() returns the first of the candidate elements that exists, so the
# script stops with a clear message if a future version renames them.
get_tab <- function(s, candidates, required = TRUE) {
  for (nm in candidates) {
    x <- s[[nm]]
    if (is.data.frame(x) && nrow(x) > 0) return(x)
  }
  if (required) stop("None of [", paste(candidates, collapse = ", "),
                     "] found in the summary. Available: ",
                     paste(names(s), collapse = ", "))
  NULL
}

# One parameter's posterior mean, 95% credible interval, ESS and R-hat
get_par <- function(tab, row) {
  if (is.null(tab) || !(row %in% rownames(tab))) {
    return(c(mean = NA, lower = NA, upper = NA, ess = NA, rhat = NA))
  }
  pick <- function(col) if (col %in% colnames(tab)) as.numeric(tab[row, col]) else NA
  c(mean = pick("Mean"), lower = pick("0.025"), upper = pick("0.975"),
    ess  = pick("ESS"),  rhat  = pick("R_hat"))
}

# One row of an inclusion table: prior probability, posterior probability and
# Bayes factor. (RoBMA's default prior inclusion probability is 0.5 for each
# component: effect, heterogeneity and publication bias.)
# When every posterior draw falls in models that include a component, the
# Bayes factor can only be bounded; RoBMA then marks it with ">" (for
# example "> 29,999"). The bound is kept so it can be displayed.
get_incl <- function(tab, row) {
  if (is.null(tab) || !(row %in% rownames(tab))) return(list(prior = NA, post = NA, bf = NA, bound = NA))
  i  <- which(rownames(tab) == row)
  bf <- tab$inclusion_BF
  bound <- attr(bf, "bound_operator")
  list(prior = if ("prior_prob" %in% names(tab)) as.numeric(tab$prior_prob[i]) else NA,
       post  = as.numeric(tab$post_prob[i]),
       bf    = as.numeric(unclass(bf))[i],
       bound = if (is.null(bound)) NA else bound[i])
}

# Formatting helpers
fmt_bf <- function(x) {
  if (is.na(x$bf)) return("--")
  val <- if (x$bf >= 100) formatC(x$bf, format = "f", digits = 0, big.mark = ",")
         else sprintf("%.2f", x$bf)
  if (!is.na(x$bound) && x$bound == ">") paste(">", val) else val
}
f3  <- function(x) if (is.na(x)) "--" else sprintf("%.3f", x)
cri <- function(lo, hi) if (is.na(lo)) "--" else sprintf("[%.3f, %.3f]", lo, hi)
f0  <- function(x) if (is.na(x)) "--" else formatC(x, format = "f", digits = 0, big.mark = ",")

# ---- Extract the tables ------------------------------------------------------
comp_rob <- get_tab(s_rob, "inclusion_components")
comp_bma <- get_tab(s_bma, "inclusion_components")
rand_rob <- get_tab(s_rob, "inclusion_random")
rand_bma <- get_tab(s_bma, "inclusion_random")
mu_rob   <- get_tab(s_rob, c("estimates_mods", "estimates"))
mu_bma   <- get_tab(s_bma, c("estimates_mods", "estimates"))
cond_rob <- get_tab(s_rob, c("estimates_mods_conditional", "estimates_conditional"), required = FALSE)
cond_bma <- get_tab(s_bma, c("estimates_mods_conditional", "estimates_conditional"), required = FALSE)
het_rob  <- get_tab(s_rob, c("estimates_random", "estimates"))
het_bma  <- get_tab(s_bma, c("estimates_random", "estimates"))
bias_rob <- get_tab(s_rob, "estimates_bias")

# ---- Sheet 1: component inclusion ---------------------------------------------
r_eff <- get_incl(comp_rob, "Effect");    b_eff <- get_incl(comp_bma, "Effect")
r_het <- get_incl(rand_rob, "tau_total"); b_het <- get_incl(rand_bma, "tau_total")
r_pb  <- get_incl(comp_rob, "Publication Bias")
comp <- data.frame(
  Component = c("Effect", "Heterogeneity", "Publication bias"),
  `RoBMA: prior probability`     = c(f3(r_eff$prior), f3(r_het$prior), f3(r_pb$prior)),
  `RoBMA: posterior probability` = c(f3(r_eff$post), f3(r_het$post), f3(r_pb$post)),
  `RoBMA: inclusion BF`          = c(fmt_bf(r_eff), fmt_bf(r_het), fmt_bf(r_pb)),
  `BMA: prior probability`       = c(f3(b_eff$prior), f3(b_het$prior), "--"),
  `BMA: posterior probability`   = c(f3(b_eff$post), f3(b_het$post), "--"),
  `BMA: inclusion BF`            = c(fmt_bf(b_eff), fmt_bf(b_het), "--"),
  check.names = FALSE, stringsAsFactors = FALSE
)

# ---- Sheet 2: estimates ------------------------------------------------------
# tau2_prop(newid) is the share of the total heterogeneity variance that lies
# between studies; the rest lies between estimates within studies.
rows <- list(
  list(label = "Mean effect, model-averaged",            tab = "mu",   par = "mu"),
  list(label = "Mean effect, conditional on an effect",  tab = "cond", par = "mu"),
  list(label = "Total heterogeneity (tau)",              tab = "het",  par = "tau_total"),
  list(label = "Share of heterogeneity between studies", tab = "het",  par = "tau2_prop(newid)")
)
est_out <- do.call(rbind, lapply(rows, function(r) {
  pr <- get_par(switch(r$tab, mu = mu_rob, cond = cond_rob, het = het_rob), r$par)
  pb <- get_par(switch(r$tab, mu = mu_bma, cond = cond_bma, het = het_bma), r$par)
  data.frame(Parameter = r$label,
             `RoBMA: mean` = f3(pr["mean"]), `RoBMA: 95% CrI` = cri(pr["lower"], pr["upper"]),
             `RoBMA: ESS` = f0(pr["ess"]),   `RoBMA: R-hat` = if (is.na(pr["rhat"])) "--" else sprintf("%.3f", pr["rhat"]),
             `BMA: mean` = f3(pb["mean"]),   `BMA: 95% CrI` = cri(pb["lower"], pb["upper"]),
             `BMA: ESS` = f0(pb["ess"]),     `BMA: R-hat` = if (is.na(pb["rhat"])) "--" else sprintf("%.3f", pb["rhat"]),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

# ---- Sheet 3: selection weights ----------------------------------------------
# omega[a,b] is the relative probability that an estimate whose one-sided
# p-value lies between a and b is reported, compared with a significant
# positive estimate (p < 0.025, fixed at 1). PET and PEESE are the regression
# coefficients of the PET and PEESE models; they are exactly zero when those
# models receive no posterior draws.
bias_rows <- rownames(bias_rob)
sel_out <- do.call(rbind, lapply(bias_rows, function(rn) {
  p <- get_par(bias_rob, rn)
  label <- if (grepl("^omega", rn)) {
    lims <- strsplit(gsub("omega\\[|\\]", "", rn), ",")[[1]]
    paste0("One-sided p from ", lims[1], " to ", lims[2])
  } else rn
  no_draws <- !grepl("^omega", rn) && (is.na(p["ess"]) || p["ess"] == 0)
  data.frame(Parameter = label,
             Mean      = if (no_draws) "0" else f3(p["mean"]),
             `95% CrI` = if (no_draws || rn == "omega[0,0.025]") "--" else cri(p["lower"], p["upper"]),
             ESS       = if (no_draws || rn == "omega[0,0.025]") "--" else f0(p["ess"]),
             `R-hat`   = if (no_draws || rn == "omega[0,0.025]") "--" else sprintf("%.3f", p["rhat"]),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
sel_out$Parameter[sel_out$Parameter == "One-sided p from 0 to 0.025"] <-
  "One-sided p from 0 to 0.025 (reference, fixed at 1)"

# ---- Sheet 4: numeric values for Protocol_Diagnostics.R ----------------------
v <- function(q, p) data.frame(Quantity = q, Estimate = p["mean"], Lower = p["lower"],
                               Upper = p["upper"], row.names = NULL, stringsAsFactors = FALSE)
values <- rbind(
  v("RoBMA mu, model-averaged",   get_par(mu_rob,   "mu")),
  v("RoBMA mu, conditional",      get_par(cond_rob, "mu")),
  v("BMA mu",                     get_par(mu_bma,   "mu")),
  v("RoBMA omega[0.025,0.05]",    get_par(bias_rob, "omega[0.025,0.05]")),
  v("RoBMA omega[0.05,0.5]",      get_par(bias_rob, "omega[0.05,0.5]")),
  v("RoBMA omega[0.5,0.95]",      get_par(bias_rob, "omega[0.5,0.95]")),
  v("RoBMA omega[0.95,0.975]",    get_par(bias_rob, "omega[0.95,0.975]")),
  v("RoBMA omega[0.975,1]",       get_par(bias_rob, "omega[0.975,1]")),
  # A Bayes factor reported only as a lower bound is stored as missing
  data.frame(Quantity = "RoBMA effect inclusion BF",
             Estimate = if (is.na(r_eff$bound)) r_eff$bf else NA,
             Lower = NA, Upper = NA, row.names = NULL)
)

write_xlsx(list(Components = comp, Estimates = est_out,
                SelectionWeights = sel_out, Values = values),
           path = here("Table_RoBMA_CHE_Intercept.xlsx"))
cat("\nSaved Table_RoBMA_CHE_Intercept.xlsx\n")
print(comp, row.names = FALSE)
print(est_out, row.names = FALSE)
print(sel_out, row.names = FALSE)

elapsed <- proc.time() - start_time
cat(sprintf("Total run time: %.1f minutes.\n", elapsed["elapsed"] / 60))

# =============================================================================
# ADAPTING THIS SCRIPT FOR A DIFFERENT EFFECT SIZE
# =============================================================================
# Fisher's z is one of six standardized effect sizes that RoBMA supports
# directly through the measure argument, with no further setup: "SMD", "ZCOR",
# "RR", "OR", "HR" and "IRR". Each has a known unit information standard
# deviation (UISD) built into the package, which RoBMA uses to scale its
# default priors.
#
# A nonstandard effect size, such as a raw regression coefficient, has no known
# UISD. It needs measure = "GEN" plus one of the following (per correspondence
# with Frantisek Bartos, RoBMA's developer, 2026-07-28):
#
#   1. If sample sizes are available for every estimate, pass them as ni:
#        RoBMA.mv(yi = yi, V = V_mat, ni = ni, measure = "GEN", ...)
#      RoBMA then estimates the UISD from ni and the sampling variances.
#
#   2. If sample sizes are available for only some estimates, compute the UISD
#      once from that subset and pass it in:
#        uisd <- estimate_unit_information_sd(sei = sei_subset, ni = ni_subset)
#        RoBMA.mv(yi = yi, V = V_mat, measure = "GEN",
#                 prior_unit_information_sd = uisd, ...)
#
# This affects only the RoBMA scripts; the metafor-based models need only yi
# and the sampling variances. The atanh()/tanh() transformations used in this
# protocol are specific to correlations and should be omitted for other
# effect sizes.
# =============================================================================
