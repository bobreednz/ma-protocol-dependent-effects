# =============================================================================
# Protocol_CHE_PETPEESE_MetaReg.R
# =============================================================================
# Purpose: the CHE-PET-PEESE meta-regression with the 18 moderators
# (Section 7.1). Replaces Protocol_Table6.R (BMS model averaging dropped).
#
# One regression answers three questions:
#   1. The corrected mean effect with the moderators held at their sample
#      means. The moderators are centered at their sample means, so the
#      intercept is the predicted effect for a study with a standard error of
#      zero and average characteristics.
#   2. Whether there is still evidence of a small-study effect once the
#      moderators are controlled for: the coefficient on the bias term.
#   3. Which study characteristics are associated with larger or smaller
#      effects: the moderator coefficients. Centering changes only the
#      intercept, not these coefficients.
#
# PET-PEESE rule (as in Section 5.2): report PEESE (bias term sez^2) if the
# PET intercept is positive with a two-sided p-value below 0.10 (a one-sided
# test at 5 percent); otherwise report PET (bias term sez).
#
# Dependence model: CHE (rho = 0.5), REML, CR2 cluster-robust standard
# errors with Satterthwaite degrees of freedom, as everywhere in the protocol.
#
# Output: Table_CHE_PETPEESE_MetaReg.xlsx
#   "Summary"    -- CHE-PET and CHE-PEESE columns: corrected mean, bias
#                   coefficient, tau, omega, observations, studies
#   "Moderators" -- moderator coefficients from the model the rule selects,
#                   with the Satterthwaite degrees of freedom of each test
#   "Values"     -- unformatted numbers, for use in later sections
# Run time: under a minute.
# =============================================================================

# -- Packages -----------------------------------------------------------------

packages <- c("here", "metafor", "clubSandwich", "writexl")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)          # resolves file paths relative to the project root (.Rproj)
library(metafor)       # rma.mv() for the CHE model
library(clubSandwich)  # coef_test() for CR2 standard errors
library(writexl)       # write_xlsx()

start_time <- proc.time()

# -- Data and covariance matrix ----------------------------------------------

DT <- readRDS(here("SCData_processed.rds"))
n_obs     <- nrow(DT)
n_studies <- length(unique(DT$newid))

rho <- 0.5   # assumed within-study correlation (Section 4.3)
V_mat <- vcalc(vi = sez^2, cluster = newid, obs = obs, rho = rho, data = DT)

# -- The 18 moderators, centered at their sample means ------------------------

moderators <- c(
  "SC1_Cognitive", "SC1_Structural", "DV_GrowthRate", "PubYear", "Published",
  "LaggedDV", "LaggedSC", "NumberSCVars", "Endog_IV", "Endog_FE",
  "CityLevel", "RegionLevel", "CountryLevel", "PanelData",
  "Reg_OECDEurope", "Reg_US", "Reg_Africa", "Reg_Asia"
)

# Labels used in the table (same wording as Table 1)
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

# Center each moderator: new variable = raw value minus its sample mean.
# as.numeric() strips any labels the variables carry from Stata.
for (m in moderators) {
  DT[[paste0(m, "_c")]] <- as.numeric(DT[[m]]) - mean(as.numeric(DT[[m]]))
}
rhs_mods <- paste(paste0(moderators, "_c"), collapse = " + ")

# -- Fit one CHE model with a given bias term --------------------------------

fit_che <- function(bias_term) {
  fml <- as.formula(paste("z ~ 1 +", bias_term, "+", rhs_mods))
  fit <- rma.mv(fml, V = V_mat, random = ~ 1 | newid/obs,
                data = DT, method = "REML")
  ct  <- as.data.frame(coef_test(fit, vcov = "CR2", cluster = DT$newid))
  # 95% confidence intervals with Satterthwaite degrees of freedom
  t_crit   <- qt(0.975, df = ct$df_Satt)
  ct$lower <- ct$beta - t_crit * ct$SE
  ct$upper <- ct$beta + t_crit * ct$SE
  list(fit = fit, ct = ct,
       tau = sqrt(fit$sigma2[1]), omega = sqrt(fit$sigma2[2]))
}

cat("Fitting CHE-PET with 18 moderators...\n")
pet   <- fit_che("sez")
cat("Fitting CHE-PEESE with 18 moderators...\n")
peese <- fit_che("I(sez^2)")

# -- Apply the PET-PEESE rule -------------------------------------------------

pet_int   <- pet$ct[pet$ct$Coef == "intrcpt", ]
use_peese <- pet_int$beta > 0 && pet_int$p_Satt < 0.10
selected  <- if (use_peese) "PEESE" else "PET"
sel       <- if (use_peese) peese else pet
cat(sprintf("PET intercept = %.3f (p = %.3f): the PET-PEESE rule selects %s.\n",
            pet_int$beta, pet_int$p_Satt, selected))

# -- Formatting helpers -------------------------------------------------------

f3    <- function(x) formatC(x, digits = 3, format = "f")
stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.10, "*", "")))
est   <- function(row) paste0(f3(row$beta), stars(row$p_Satt))
se    <- function(row) paste0("(", f3(row$SE), ")")
ci    <- function(row) paste0("[", f3(row$lower), ", ", f3(row$upper), "]")
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", f3(p))

# -- Sheet 1: Summary (both bias terms) ---------------------------------------

summary_column <- function(m, bias_name) {
  int  <- m$ct[m$ct$Coef == "intrcpt", ]
  bias <- m$ct[m$ct$Coef == bias_name, ]
  c(est(int), se(int), ci(int), est(bias), se(bias),
    f3(m$tau), f3(m$omega), n_obs, n_studies)
}

summary_tbl <- data.frame(
  Parameter = c("Mean effect", "(SE)", "95% CI",
                "Bias coefficient", "(SE)",
                "tau", "omega", "Observations", "Studies"),
  `CHE-PET`   = summary_column(pet,   "sez"),
  `CHE-PEESE` = summary_column(peese, "I(sez^2)"),
  check.names = FALSE, stringsAsFactors = FALSE
)

# -- Sheet 2: Moderator coefficients (selected model) -------------------------

mods_ct <- sel$ct[match(paste0(moderators, "_c"), sel$ct$Coef), ]
moderator_tbl <- data.frame(
  Moderator   = unname(mod_labels[moderators]),
  Coefficient = est(mods_ct),
  `(SE)`      = se(mods_ct),
  # Satterthwaite df: tests with fewer than about 4 are unreliable (Tipton, 2015)
  df          = formatC(mods_ct$df_Satt, digits = 1, format = "f"),
  `p-value`   = fmt_p(mods_ct$p_Satt),
  check.names = FALSE, stringsAsFactors = FALSE
)

# -- Sheet 3: unformatted values ---------------------------------------------

value_row <- function(model, quantity, row) {
  data.frame(Model = model, Quantity = quantity, Estimate = row$beta,
             SE = row$SE, Lower = row$lower, Upper = row$upper, p = row$p_Satt)
}
values <- rbind(
  value_row("CHE-PET",   "Mean effect", pet$ct[pet$ct$Coef == "intrcpt", ]),
  value_row("CHE-PET",   "Bias coefficient", pet$ct[pet$ct$Coef == "sez", ]),
  value_row("CHE-PEESE", "Mean effect", peese$ct[peese$ct$Coef == "intrcpt", ]),
  value_row("CHE-PEESE", "Bias coefficient", peese$ct[peese$ct$Coef == "I(sez^2)", ]),
  value_row(paste0("CHE-", selected), unname(mod_labels[moderators]), mods_ct)
)
values$Selected <- values$Model == paste0("CHE-", selected)

# -- Save ---------------------------------------------------------------------

write_xlsx(list(Summary = summary_tbl, Moderators = moderator_tbl, Values = values),
           path = here("Table_CHE_PETPEESE_MetaReg.xlsx"))

print(summary_tbl, row.names = FALSE)
print(moderator_tbl, row.names = FALSE)
cat(sprintf("\nSaved Table_CHE_PETPEESE_MetaReg.xlsx. Run time: %.1f seconds.\n",
            (proc.time() - start_time)["elapsed"]))
