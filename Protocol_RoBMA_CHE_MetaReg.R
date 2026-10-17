# =============================================================================
# Protocol_RoBMA_CHE_MetaReg.R
# =============================================================================
# Purpose: RoBMA with the 18 moderators, and its best-practice predictions
# (Section 7.3).
#
#   RoBMA.mv: the CHE RoBMA of Section 5.4 with the 18 moderators in the mean.
#             RoBMA averages over models with and without each moderator, so
#             every moderator gets its own inclusion Bayes factor (prior
#             inclusion probability 0.5), as the effect, heterogeneity and
#             publication bias do in Section 5.4.
#   BMA.mv:   the same model without the publication-bias component (as in
#             Step 1 of Section 5.4). Comparing its inclusion Bayes factors
#             with RoBMA's shows which moderators the data support before and
#             after the correction.
#
# What the script reports:
#   1. The corrected mean at average characteristics: pooled_effect() predicts
#      the mean effect for every estimate from its own moderator values and
#      averages the predictions. It plays the role of the centered intercept
#      in Sections 7.1 and 7.2.
#   2. Whether selection is still present with the moderators included: the
#      inclusion Bayes factor for publication bias and the selection weights.
#   3. Which moderators belong in the model: the inclusion Bayes factors and
#      the model-averaged coefficients.
#   4. Best-practice predictions from predict(): the mean effect at the
#      best-practice values (type = "terms", with a credible interval) and
#      the distribution of true effects there (type = "estimate", giving the
#      prediction interval; it does not include sampling error).
#
# Coding of the moderators (as in the saved fits):
#   - Binary moderators enter as factors with treatment contrasts, so each
#     coefficient is the difference between estimates with and without the
#     characteristic.
#   - RoBMA standardizes the continuous moderators (PubYear, NumberSCVars)
#     internally. summary() reports their coefficients in the original units
#     (per year, per additional variable), which is what the table shows.
#   - The intercept is NOT reported. In the original units it refers to
#     PubYear = 0; even with standardized coefficients it is the effect for a
#     study at the reference level of every binary moderator, which is not a
#     meaningful profile here. RoBMA's "Effect" Bayes factor refers to this
#     intercept, so it is not reported either.
#
# SOFTWARE: RoBMA 4.1 (development version; see Section 5.4) and JAGS.
#
# SAVED FITS: refitting takes about 5 hours for RoBMA.mv and about 30 minutes
# for BMA.mv. By default the script loads the saved fits (refit <- FALSE).
#
# Outputs:
#   fit_RoBMA_CHE_reg18.rds, fit_BMA_CHE_reg18.rds   fitted models
#   Table_RoBMA_CHE_MetaReg.xlsx
#     "Summary"               mean at average characteristics, heterogeneity,
#                             publication-bias Bayes factor (RoBMA and BMA)
#     "Moderators"            inclusion Bayes factors and model-averaged
#                             coefficients (RoBMA), inclusion Bayes factors (BMA)
#     "SelectionWeights"      RoBMA selection weights and PET/PEESE
#     "BP1 - Non-OECD Europe" best-practice prediction, Fisher's z and PCC
#     "BP2 - OECD Europe"     the same for OECD/Europe
#     "Notes"                 the best-practice values, one note per scenario
#     "Convergence"           ESS and R-hat for every reported parameter
#     "Values"                unformatted numbers, for use in Section 8
# Run time with the saved fits: a few minutes (mostly predict()).
# =============================================================================

start_time <- proc.time()   # record start time so run time can be reported

refit <- FALSE              # TRUE re-estimates both models (about 5.5 hours)

# ---- Packages ----------------------------------------------------------------
library(here)       # resolves file paths relative to the project root (.Rproj)
library(metafor)    # vcalc() for the sampling covariance matrix
library(RoBMA)      # RoBMA.mv(), BMA.mv(), pooled_effect(); version 4.1 or later
library(posterior)  # as_draws(), to extract posterior draws
library(writexl)    # write_xlsx()

# The CHE models need RoBMA 4.1 or later: RoBMA.mv() is not in the CRAN version (4.0.0).
if (packageVersion("RoBMA") < "4.1") {
  stop("This script needs the development version of RoBMA (4.1 or later). ",
       "Run Protocol_Install.R, or see Section 5.4 of the protocol.")
}

cat("RoBMA version:", as.character(packageVersion("RoBMA")), "\n")

rob_file <- here("fit_RoBMA_CHE_reg18.rds")
bma_file <- here("fit_BMA_CHE_reg18.rds")

# ---- Data ----------------------------------------------------------------------
dat <- readRDS(here("SCData_processed.rds"))
dat$yi <- atanh(dat$pcc)   # Fisher's z

# The 0/1 moderators must be factors: RoBMA treats numeric predictors as
# continuous and standardizes them. PubYear and NumberSCVars stay numeric.
binary_mods <- c("SC1_Cognitive", "SC1_Structural", "DV_GrowthRate", "Published",
                 "LaggedDV", "LaggedSC", "Endog_IV", "Endog_FE", "CityLevel",
                 "RegionLevel", "CountryLevel", "PanelData", "Reg_OECDEurope",
                 "Reg_US", "Reg_Africa", "Reg_Asia")
for (v in binary_mods) dat[[v]] <- as.factor(dat[[v]])

# Same 18 moderators, in the same order, as in Sections 7.1 and 7.2
moderators <- c(
  "SC1_Cognitive", "SC1_Structural", "DV_GrowthRate", "PubYear", "Published",
  "LaggedDV", "LaggedSC", "NumberSCVars", "Endog_IV", "Endog_FE",
  "CityLevel", "RegionLevel", "CountryLevel", "PanelData",
  "Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia"
)
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
continuous <- c("PubYear", "NumberSCVars")

# ---- Fit (only if requested or if the saved fits are missing) -------------
if (refit || !file.exists(rob_file) || !file.exists(bma_file)) {

  # Correlated sampling errors (the "C" in CHE), as in the frequentist models
  V_mat <- vcalc(vi = sez^2, cluster = newid, obs = obs, rho = 0.5, data = dat)
  mod_formula <- as.formula(paste("~", paste(moderators, collapse = " + ")))

  cat("Fitting the 18-moderator CHE RoBMA (about 5 hours)...\n")
  fit_rob <- RoBMA.mv(
    yi = yi, V = V_mat, measure = "ZCOR",
    mods = mod_formula,
    set_contrast_factor_predictors = "treatment",   # dummy coding, as in 8.1 and 8.2
    random = ~ 1 | newid/obs,                       # study and estimate random effects
    sample = 5000, burnin = 5000, adapt = 5000, thin = 3,
    parallel = TRUE, seed = 1, data = dat)
  saveRDS(fit_rob, rob_file)

  cat("Fitting the 18-moderator CHE BMA, no bias component (about 30 minutes)...\n")
  fit_bma <- BMA.mv(
    yi = yi, V = V_mat, measure = "ZCOR",
    mods = mod_formula,
    set_contrast_factor_predictors = "treatment",
    random = ~ 1 | newid/obs,
    # Longer chains (thin = 15): with thin = 3, three parameters fell short
    # of an ESS of 500
    sample = 5000, burnin = 5000, adapt = 5000, thin = 15,
    parallel = TRUE, seed = 1, data = dat)
  saveRDS(fit_bma, bma_file)

} else {
  cat("Loading saved fits (set refit <- TRUE to re-estimate)...\n")
  fit_rob <- readRDS(rob_file)
  fit_bma <- readRDS(bma_file)
}

# ---- Summaries ---------------------------------------------------------------
# summary() reports model-averaged estimates, with continuous moderators in
# their original units. Printing it also shows R-hat and ESS for every
# parameter.
s_rob <- summary(fit_rob)
s_bma <- summary(fit_bma)
print(s_rob)
print(s_bma)
cat("\nElements of the RoBMA summary object:", paste(names(s_rob), collapse = ", "), "\n")

# ---- Helpers for reading the summary tables ----------------------------------
# The summary object is a list of tables. Rather than rely on the names of the
# list elements, find_tab() returns the first table that has a given row and
# a given column. This keeps the script working if a later version of RoBMA
# renames the elements.
find_tab <- function(s, row, col) {
  for (nm in names(s)) {
    x <- s[[nm]]
    if (is.data.frame(x) && row %in% rownames(x) && col %in% colnames(x)) return(x)
  }
  stop("No table in the summary has row '", row, "' and column '", col,
       "'. Elements: ", paste(names(s), collapse = ", "))
}

# Row names of the coefficients: factors carry a level suffix, e.g. "Endog_FE[1]"
coef_row <- function(tab, v) {
  hit <- rownames(tab)[rownames(tab) == v | startsWith(rownames(tab), paste0(v, "["))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# One parameter's posterior mean, 95% credible interval, ESS and R-hat
get_par <- function(tab, row) {
  if (is.na(row) || !(row %in% rownames(tab))) {
    return(c(mean = NA, lower = NA, upper = NA, ess = NA, rhat = NA))
  }
  pick <- function(col) if (col %in% colnames(tab)) as.numeric(tab[row, col]) else NA
  c(mean = pick("Mean"), lower = pick("0.025"), upper = pick("0.975"),
    ess  = pick("ESS"),  rhat  = pick("R_hat"))
}

# Inclusion Bayes factor. When every posterior draw falls in models with the
# component, RoBMA can give only a lower bound and marks it with ">".
get_bf <- function(tab, row) {
  if (!(row %in% rownames(tab))) return(list(bf = NA, bound = NA))
  i  <- which(rownames(tab) == row)
  bf <- tab$inclusion_BF
  bound <- attr(bf, "bound_operator")
  list(bf = as.numeric(unclass(bf))[i], bound = if (is.null(bound)) NA else bound[i])
}

# Posterior mean and 95% interval from draws or, failing that, from summary()
draw_summary <- function(obj) {
  x <- tryCatch(as.numeric(as_draws(obj)), error = function(e) NULL)
  if (!is.null(x) && length(x) > 1) {
    return(c(mean = mean(x), lower = unname(quantile(x, 0.025)),
             upper = unname(quantile(x, 0.975))))
  }
  s <- as.data.frame(summary(obj))
  c(mean = as.numeric(s[["Mean"]][1]), lower = as.numeric(s[["0.025"]][1]),
    upper = as.numeric(s[["0.975"]][1]))
}

# Formatting helpers
fmt_bf <- function(x) {
  if (is.na(x$bf)) return("--")
  val <- if (x$bf >= 100) formatC(x$bf, format = "f", digits = 0, big.mark = ",")
         else sprintf("%.2f", x$bf)
  if (!is.na(x$bound) && x$bound == ">") paste(">", val) else val
}
# no_neg_zero() turns "-0.000" (a tiny negative number) into "0.000"
no_neg_zero <- function(txt) sub("^-(0\\.0+)$", "\\1", txt)
f3  <- function(x) ifelse(is.na(x), "--", no_neg_zero(sprintf("%.3f", x)))
f4  <- function(x) ifelse(is.na(x), "--", no_neg_zero(sprintf("%.4f", x)))
cri <- function(lo, hi, f = f3) ifelse(is.na(lo), "--", paste0("[", f(lo), ", ", f(hi), "]"))
est_cri <- function(p) if (is.na(p["mean"])) "--" else paste(f3(p["mean"]), cri(p["lower"], p["upper"]))

# ---- Extract the tables ------------------------------------------------------
incl_mods_rob <- find_tab(s_rob, "NumberSCVars", "inclusion_BF")   # moderator inclusion
incl_mods_bma <- find_tab(s_bma, "NumberSCVars", "inclusion_BF")
coef_rob      <- find_tab(s_rob, "PubYear", "Mean")                # moderator coefficients
coef_bma      <- find_tab(s_bma, "PubYear", "Mean")
het_rob       <- find_tab(s_rob, "tau_total", "Mean")
het_bma       <- find_tab(s_bma, "tau_total", "Mean")
comp_rob      <- find_tab(s_rob, "Publication Bias", "inclusion_BF")
bias_rob      <- find_tab(s_rob, "omega[0.5,0.95]", "Mean")

# ---- Pooled effect: the mean at average characteristics ----------------------
cat("\nComputing pooled effects...\n")
pe_rob <- draw_summary(pooled_effect(fit_rob))
pe_bma <- draw_summary(pooled_effect(fit_bma))

# ---- Sheet 1: Summary --------------------------------------------------------
# tau2_prop(newid) is the share of the total heterogeneity variance that lies
# between studies; the rest lies between estimates within studies.
tau_rob <- get_par(het_rob, "tau_total");        tau_bma <- get_par(het_bma, "tau_total")
shr_rob <- get_par(het_rob, "tau2_prop(newid)"); shr_bma <- get_par(het_bma, "tau2_prop(newid)")
pb_rob  <- get_bf(comp_rob, "Publication Bias")

summary_tbl <- data.frame(
  Parameter = c("Mean effect at average characteristics",
                "Total heterogeneity (tau)",
                "Share of heterogeneity between studies",
                "Publication bias, inclusion BF"),
  RoBMA = c(est_cri(pe_rob), est_cri(tau_rob), f3(shr_rob["mean"]), fmt_bf(pb_rob)),
  `BMA (no correction)` = c(est_cri(pe_bma), est_cri(tau_bma), f3(shr_bma["mean"]), "--"),
  check.names = FALSE, stringsAsFactors = FALSE
)

# ---- Sheet 2: Moderators -----------------------------------------------------
# The coefficients are model-averaged: they average over models with and
# without the moderator (where the coefficient is zero), so a moderator with
# little support has a coefficient close to zero, and a credible interval can
# end at exactly zero.
mod_rows <- lapply(moderators, function(v) {
  p  <- get_par(coef_rob, coef_row(coef_rob, v))
  fx <- if (v %in% continuous) f4 else f3          # per-unit effects are small
  data.frame(
    Moderator                  = unname(mod_labels[v]),
    `RoBMA: coefficient`       = fx(p["mean"]),
    `RoBMA: 95% CrI`           = cri(p["lower"], p["upper"], fx),
    # The two Bayes factors sit side by side, so they can be compared directly
    `RoBMA: inclusion BF`      = fmt_bf(get_bf(incl_mods_rob, v)),
    `BMA: inclusion BF`        = fmt_bf(get_bf(incl_mods_bma, v)),
    check.names = FALSE, stringsAsFactors = FALSE)
})
moderator_tbl <- do.call(rbind, mod_rows)

# ---- Sheet 3: Selection weights (as in Section 5.4, Panel C) -----------------
sel_tbl <- do.call(rbind, lapply(rownames(bias_rob), function(rn) {
  p <- get_par(bias_rob, rn)
  label <- if (grepl("^omega", rn)) {
    lims <- strsplit(gsub("omega\\[|\\]", "", rn), ",")[[1]]
    paste0("One-sided p from ", lims[1], " to ", lims[2])
  } else rn
  ref      <- rn == "omega[0,0.025]"
  no_draws <- !grepl("^omega", rn) && (is.na(p["ess"]) || p["ess"] == 0)
  data.frame(Parameter        = if (ref) paste(label, "(reference, fixed at 1)") else label,
             `Posterior mean` = if (no_draws) "0" else f3(p["mean"]),
             `95% CrI`        = if (no_draws || ref) "--" else cri(p["lower"], p["upper"]),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

# ---- Sheets 4 and 5: best-practice predictions -------------------------------
# Same best-practice values as Protocol_CHE_PETPEESE_BestPractice.R. BP1 = non-OECD/Europe (all
# region dummies 0); BP2 = OECD/Europe. Factor levels are copied from the data
# so that predict() codes the new data exactly as the fitted model does.
make_bp <- function(oecd) {
  vals <- list(SC1_Cognitive = 1, SC1_Structural = 0, DV_GrowthRate = 1,
               Published = 1, LaggedDV = 0, LaggedSC = 0, Endog_IV = 1,
               Endog_FE = 1, CityLevel = 0, RegionLevel = 0, CountryLevel = 1,
               PanelData = 1, Reg_OECDEurope = oecd, Reg_US = 0,
               Reg_Africa = 0, Reg_Asia = 0)
  nd <- data.frame(PubYear = mean(dat$PubYear), NumberSCVars = 1)
  for (v in names(vals)) nd[[v]] <- factor(vals[[v]], levels = levels(dat[[v]]))
  nd
}

bp_panel <- function(oecd) {
  nd <- make_bp(oecd)
  # predict(type = "estimate") simulates the random effects, so fix the seed
  # to make the prediction interval reproducible
  set.seed(20261017)
  # type = "terms": the mean effect at these values (credible interval)
  m  <- draw_summary(predict(fit_rob, newdata = nd, type = "terms",    bias_adjusted = TRUE))
  # type = "estimate": a true effect at these values (prediction interval);
  # adds the heterogeneity but not the sampling error of a new estimate
  pr <- draw_summary(predict(fit_rob, newdata = nd, type = "estimate", bias_adjusted = TRUE))
  list(
    display = data.frame(
      `Effect Size`     = c("Fisher's z", "PCC"),
      `Mean Prediction` = c(f3(m["mean"]), f3(tanh(m["mean"]))),
      `95% CrI`         = c(cri(m["lower"], m["upper"]), cri(tanh(m["lower"]), tanh(m["upper"]))),
      `95% PI`          = c(cri(pr["lower"], pr["upper"]), cri(tanh(pr["lower"]), tanh(pr["upper"]))),
      check.names = FALSE, stringsAsFactors = FALSE),
    values = c(m, pi_lo = unname(pr["lower"]), pi_hi = unname(pr["upper"]))
  )
}
cat("Computing best-practice predictions...\n")
bp1 <- bp_panel(0)
bp2 <- bp_panel(1)

bp_note <- function(oecd) {
  paste0(
    "\"Best Practice\" predictions are derived from the RoBMA meta-regression with 18 moderators (Section 7.3). ",
    "Predictions are evaluated at SC1_Cognitive = 1, SC1_Structural = 0, ",
    "SC1_Other = 0 (reference), DV_GrowthRate = 1, DV_GDPLevel = 0 (reference), PubYear = ",
    round(mean(dat$PubYear), 1), " (sample mean), Published = 1, LaggedDV = 0, ",
    "LaggedSC = 0, NumberSCVars = 1, Endog_IV = 1, Endog_FE = 1, CityLevel = 0, RegionLevel = 0, ",
    "CountryLevel = 1, PanelData = 1, Reg_OECDEurope = ", oecd,
    ", Reg_US = 0, Reg_Africa = 0, Reg_Asia = 0. ",
    "The prediction is model-averaged and corrected for publication bias. CrI = credible interval. ",
    "The prediction interval is the range of true effects at these values; it does not include sampling error."
  )
}
notes_sheet <- data.frame(Scenario = c("BP1 - Non-OECD/Europe", "BP2 - OECD/Europe"),
                          Note = c(bp_note(0), bp_note(1)), stringsAsFactors = FALSE)

# ---- Sheet 6: Convergence ----------------------------------------------------
# ESS and R-hat for every reported parameter. RoBMA's default criteria are an
# R-hat below 1.05 and an ESS of at least 500.
conv_rows <- function(model, coef_tab, het_tab, bias_tab = NULL) {
  pars <- c(lapply(moderators, function(v) list(label = unname(mod_labels[v]), tab = coef_tab,
                                                row = coef_row(coef_tab, v))),
            list(list(label = "Total heterogeneity (tau)", tab = het_tab, row = "tau_total")))
  if (!is.null(bias_tab)) {
    pars <- c(pars, lapply(setdiff(grep("^omega", rownames(bias_tab), value = TRUE), "omega[0,0.025]"),
                           function(rn) list(label = rn, tab = bias_tab, row = rn)))
  }
  do.call(rbind, lapply(pars, function(p) {
    g <- get_par(p$tab, p$row)
    data.frame(Model = model, Parameter = p$label, ESS = unname(g["ess"]),
               R_hat = unname(g["rhat"]), stringsAsFactors = FALSE)
  }))
}
convergence <- rbind(conv_rows("RoBMA", coef_rob, het_rob, bias_rob),
                     conv_rows("BMA", coef_bma, het_bma))
convergence$Meets_criteria <- !is.na(convergence$ESS) & convergence$ESS >= 500 &
                              !is.na(convergence$R_hat) & convergence$R_hat < 1.05

# ---- Sheet 7: unformatted values for Section 8 --------------------------------
v <- function(q, p, pi = c(NA, NA)) data.frame(
  Quantity = q, Estimate = unname(p["mean"]), Lower = unname(p["lower"]),
  Upper = unname(p["upper"]), PI_lower = pi[1], PI_upper = pi[2], stringsAsFactors = FALSE)
values <- rbind(
  v("RoBMA mean at average characteristics", pe_rob),
  v("BMA mean at average characteristics",   pe_bma),
  v("RoBMA BP1 (non-OECD/Europe)", bp1$values, bp1$values[c("pi_lo", "pi_hi")]),
  v("RoBMA BP2 (OECD/Europe)",     bp2$values, bp2$values[c("pi_lo", "pi_hi")])
)

# ---- Save ---------------------------------------------------------------------
write_xlsx(list(Summary                 = summary_tbl,
                Moderators              = moderator_tbl,
                SelectionWeights        = sel_tbl,
                "BP1 - Non-OECD Europe" = bp1$display,
                "BP2 - OECD Europe"     = bp2$display,
                Notes                   = notes_sheet,
                Convergence             = convergence,
                Values                  = values),
           path = here("Table_RoBMA_CHE_MetaReg.xlsx"))

print(summary_tbl, row.names = FALSE)
print(moderator_tbl, row.names = FALSE)
print(sel_tbl, row.names = FALSE)
cat("\nBP1 (non-OECD/Europe):\n"); print(bp1$display, row.names = FALSE)
cat("\nBP2 (OECD/Europe):\n");     print(bp2$display, row.names = FALSE)
cat("\nParameters below the convergence criteria:\n")
print(convergence[!convergence$Meets_criteria, ], row.names = FALSE)
cat(sprintf("\nSaved Table_RoBMA_CHE_MetaReg.xlsx. Run time: %.1f minutes.\n",
            (proc.time() - start_time)["elapsed"] / 60))
