setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("metafor", "clubSandwich", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(metafor)       # meta-analytic models via rma.mv()
library(clubSandwich)  # CR2 clustered standard errors
library(tidyverse)     # data manipulation
library(writexl)       # export to Excel

# ============================================================
# Load data
# ============================================================

DT <- readRDS("SCData_processed.rds")

n_obs     <- nrow(DT)
n_studies <- length(unique(DT$newid))

# ============================================================
# Impute block-diagonal covariance matrix (for CE and CHE)
# ============================================================

rho <- 0.5

V_mat <- vcalc(
  vi      = sez^2,
  cluster = newid,
  obs     = obs,
  rho     = rho,
  data    = DT
)

# ============================================================
# Center all control variables at their sample means.
#
# This is the key step for Table 5. With centered controls,
# the intercept in each model directly gives the "Effect
# beyond bias": the estimated mean effect when sez = 0 (no
# publication bias) and all controls are at their means
# (centered values = 0). This means we can read the "Effect
# beyond bias" and its CR2 standard error directly from
# coef_test(), with no post-estimation prediction required.
# ============================================================

DT <- DT %>%
  mutate(
    SC1_Cognitive_c  = SC1_Cognitive  - mean(SC1_Cognitive),
    SC1_Structural_c = SC1_Structural - mean(SC1_Structural),
    DV_GrowthRate_c  = DV_GrowthRate  - mean(DV_GrowthRate),
    PubYear_c        = PubYear        - mean(PubYear),
    Published_c      = Published      - mean(Published),
    LaggedDV_c       = LaggedDV       - mean(LaggedDV),
    LaggedSC_c       = LaggedSC       - mean(LaggedSC),
    NumberSCVars_c   = NumberSCVars   - mean(NumberSCVars),
    Endog_IV_c       = Endog_IV       - mean(Endog_IV),
    Endog_FE_c       = Endog_FE       - mean(Endog_FE),
    CityLevel_c      = CityLevel      - mean(CityLevel),
    RegionLevel_c    = RegionLevel    - mean(RegionLevel),
    CountryLevel_c   = CountryLevel   - mean(CountryLevel),
    PanelData_c      = PanelData      - mean(PanelData),
    Reg_OECDEurope_c = Reg_OECDEurope - mean(Reg_OECDEurope),
    Reg_US_c         = Reg_US         - mean(Reg_US),
    Reg_Africa_c     = Reg_Africa     - mean(Reg_Africa),
    Reg_Asia_c       = Reg_Asia       - mean(Reg_Asia)
  )

# ============================================================
# Helper functions
# ============================================================

# Significance stars
stars <- function(p) {
  case_when(p < 0.01 ~ "***", p < 0.05 ~ "**", p < 0.10 ~ "*", TRUE ~ "")
}

# Format LR test p-value
fmt_p <- function(p) {
  if (p < 0.0001) "p < .0001" else sprintf("p = %.3f", p)
}

# Format variance component as SD (square root of sigma2)
fmt_vc <- function(sigma2) formatC(sqrt(sigma2), digits = 3, format = "f")

# Format estimate with significance stars
fmt_est <- function(est, p) paste0(formatC(est, digits = 3, format = "f"), stars(p))

# Format standard error in parentheses
fmt_se <- function(se) paste0("(", formatC(se, digits = 3, format = "f"), ")")

# Extract a named row from coef_test() output
get_row <- function(ct, coef_name) ct[ct$Coef == coef_name, ]

# Format a 95% confidence interval for the intercept ("Effect
# beyond bias"), using the Satterthwaite df from coef_test() --
# the same construction used in Table 4 and Table 7, so every
# confidence interval in the protocol is built the same way.
fmt_ci <- function(row) {
  t_crit <- qt(0.975, df = row$df_Satt)
  lo <- row$beta - t_crit * row$SE
  hi <- row$beta + t_crit * row$SE
  paste0("[", formatC(lo, digits = 3, format = "f"), ", ",
              formatC(hi, digits = 3, format = "f"), "]")
}

# ============================================================
# PANEL A: sez only (no control variables)
#
# "Bias"              = coefficient on sez
# "Effect beyond bias" = intercept (the constant when sez = 0)
# ============================================================

# --- REML models ---

pA_FE  <- rma.mv(z ~ 1 + sez,
                  V = sez^2,
                  data = DT, method = "FE")

pA_RE  <- rma.mv(z ~ 1 + sez,
                  V = sez^2,
                  random = ~ 1 | obs,
                  data = DT, method = "REML")

pA_CE  <- rma.mv(z ~ 1 + sez,
                  V = V_mat,
                  random = ~ 1 | newid,
                  data = DT, method = "REML")

pA_HE  <- rma.mv(z ~ 1 + sez,
                  V = sez^2,
                  random = ~ 1 | newid/obs,
                  data = DT, method = "REML")

pA_CHE <- rma.mv(z ~ 1 + sez,
                  V = V_mat,
                  random = ~ 1 | newid/obs,
                  data = DT, method = "REML")

# --- ML models for LR tests ---

pA_RE_ML  <- rma.mv(z ~ 1 + sez, V = sez^2, random = ~ 1 | obs,       data = DT, method = "ML")
pA_CE_ML  <- rma.mv(z ~ 1 + sez, V = V_mat, random = ~ 1 | newid,     data = DT, method = "ML")
pA_HE_ML  <- rma.mv(z ~ 1 + sez, V = sez^2, random = ~ 1 | newid/obs, data = DT, method = "ML")
pA_CHE_ML <- rma.mv(z ~ 1 + sez, V = V_mat, random = ~ 1 | newid/obs, data = DT, method = "ML")

# --- CR2 standard errors ---

ct_pA_FE  <- coef_test(pA_FE,  cluster = DT$newid, vcov = "CR2")
ct_pA_RE  <- coef_test(pA_RE,  cluster = DT$newid, vcov = "CR2")
ct_pA_CE  <- coef_test(pA_CE,  cluster = DT$newid, vcov = "CR2")
ct_pA_HE  <- coef_test(pA_HE,  cluster = DT$newid, vcov = "CR2")
ct_pA_CHE <- coef_test(pA_CHE, cluster = DT$newid, vcov = "CR2")

# --- LR tests ---

lrt_pA_RE  <- anova(pA_FE,    pA_RE_ML)   # RE vs FE
lrt_pA_HE  <- anova(pA_RE_ML, pA_HE_ML)   # HE vs RE
lrt_pA_CHE <- anova(pA_CE_ML, pA_CHE_ML)  # CHE vs CE

# --- Assemble Panel A table ---

panel_A <- tibble(
  Variable = c("Effect beyond bias", "(SE)", "95% CI", "Bias", "(SE)",
               "Observations", "Studies", "rho (assumed)",
               "tau", "omega", "LR test"),

  # rho does not apply to FE, RE, or HE (no imputed V matrix)
  FE = c(
    fmt_est(get_row(ct_pA_FE, "intrcpt")$beta,  get_row(ct_pA_FE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pA_FE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pA_FE, "intrcpt")),
    fmt_est(get_row(ct_pA_FE, "sez")$beta,     get_row(ct_pA_FE, "sez")$p_Satt),
    fmt_se( get_row(ct_pA_FE, "sez")$SE),
    n_obs, n_studies, "--", "--", "--", "--"
  ),

  RE = c(
    fmt_est(get_row(ct_pA_RE, "intrcpt")$beta,  get_row(ct_pA_RE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pA_RE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pA_RE, "intrcpt")),
    fmt_est(get_row(ct_pA_RE, "sez")$beta,     get_row(ct_pA_RE, "sez")$p_Satt),
    fmt_se( get_row(ct_pA_RE, "sez")$SE),
    n_obs, n_studies, "--",
    fmt_vc(pA_RE$sigma2), "--",
    fmt_p(lrt_pA_RE$pval)
  ),

  # CE and CHE use the imputed V matrix, so rho is fixed at the
  # value assumed above, not estimated.
  CE = c(
    fmt_est(get_row(ct_pA_CE, "intrcpt")$beta,  get_row(ct_pA_CE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pA_CE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pA_CE, "intrcpt")),
    fmt_est(get_row(ct_pA_CE, "sez")$beta,     get_row(ct_pA_CE, "sez")$p_Satt),
    fmt_se( get_row(ct_pA_CE, "sez")$SE),
    n_obs, n_studies,
    formatC(rho, digits = 3, format = "f"),
    fmt_vc(pA_CE$sigma2), "--", "--"
  ),

  HE = c(
    fmt_est(get_row(ct_pA_HE, "intrcpt")$beta,  get_row(ct_pA_HE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pA_HE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pA_HE, "intrcpt")),
    fmt_est(get_row(ct_pA_HE, "sez")$beta,     get_row(ct_pA_HE, "sez")$p_Satt),
    fmt_se( get_row(ct_pA_HE, "sez")$SE),
    n_obs, n_studies, "--",
    fmt_vc(pA_HE$sigma2[1]), fmt_vc(pA_HE$sigma2[2]),
    fmt_p(lrt_pA_HE$pval)
  ),

  CHE = c(
    fmt_est(get_row(ct_pA_CHE, "intrcpt")$beta,  get_row(ct_pA_CHE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pA_CHE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pA_CHE, "intrcpt")),
    fmt_est(get_row(ct_pA_CHE, "sez")$beta,     get_row(ct_pA_CHE, "sez")$p_Satt),
    fmt_se( get_row(ct_pA_CHE, "sez")$SE),
    n_obs, n_studies,
    formatC(rho, digits = 3, format = "f"),
    fmt_vc(pA_CHE$sigma2[1]), fmt_vc(pA_CHE$sigma2[2]),
    fmt_p(lrt_pA_CHE$pval)
  )
)

# ============================================================
# PANEL B: sez + full set of centered control variables
#
# Because controls are centered, the intercept is again the
# "Effect beyond bias": estimated mean effect when sez = 0
# and all controls equal their sample means.
# ============================================================

formula_B <- z ~ 1 + sez +
  SC1_Cognitive_c + SC1_Structural_c + DV_GrowthRate_c +
  PubYear_c + Published_c + LaggedDV_c + LaggedSC_c + NumberSCVars_c +
  Endog_IV_c + Endog_FE_c + CityLevel_c + RegionLevel_c + CountryLevel_c +
  PanelData_c + Reg_OECDEurope_c + Reg_US_c + Reg_Africa_c + Reg_Asia_c

# --- REML models ---

pB_FE  <- rma.mv(formula_B, V = sez^2, data = DT, method = "FE")

pB_RE  <- rma.mv(formula_B,
                  V = sez^2,
                  random = ~ 1 | obs,
                  data = DT, method = "REML")

pB_CE  <- rma.mv(formula_B,
                  V = V_mat,
                  random = ~ 1 | newid,
                  data = DT, method = "REML")

pB_HE  <- rma.mv(formula_B,
                  V = sez^2,
                  random = ~ 1 | newid/obs,
                  data = DT, method = "REML")

pB_CHE <- rma.mv(formula_B,
                  V = V_mat,
                  random = ~ 1 | newid/obs,
                  data = DT, method = "REML")

# --- ML models for LR tests ---

pB_RE_ML  <- rma.mv(formula_B, V = sez^2, random = ~ 1 | obs,       data = DT, method = "ML")
pB_CE_ML  <- rma.mv(formula_B, V = V_mat, random = ~ 1 | newid,     data = DT, method = "ML")
pB_HE_ML  <- rma.mv(formula_B, V = sez^2, random = ~ 1 | newid/obs, data = DT, method = "ML")
pB_CHE_ML <- rma.mv(formula_B, V = V_mat, random = ~ 1 | newid/obs, data = DT, method = "ML")

# --- CR2 standard errors ---

ct_pB_FE  <- coef_test(pB_FE,  cluster = DT$newid, vcov = "CR2")
ct_pB_RE  <- coef_test(pB_RE,  cluster = DT$newid, vcov = "CR2")
ct_pB_CE  <- coef_test(pB_CE,  cluster = DT$newid, vcov = "CR2")
ct_pB_HE  <- coef_test(pB_HE,  cluster = DT$newid, vcov = "CR2")
ct_pB_CHE <- coef_test(pB_CHE, cluster = DT$newid, vcov = "CR2")

# --- LR tests ---

lrt_pB_RE  <- anova(pB_FE,    pB_RE_ML)   # RE vs FE
lrt_pB_HE  <- anova(pB_RE_ML, pB_HE_ML)   # HE vs RE
lrt_pB_CHE <- anova(pB_CE_ML, pB_CHE_ML)  # CHE vs CE

# --- Assemble Panel B table ---

panel_B <- tibble(
  Variable = c("Effect beyond bias", "(SE)", "95% CI", "Bias", "(SE)",
               "Observations", "Studies", "rho (assumed)",
               "tau", "omega", "LR test"),

  FE = c(
    fmt_est(get_row(ct_pB_FE, "intrcpt")$beta,  get_row(ct_pB_FE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pB_FE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pB_FE, "intrcpt")),
    fmt_est(get_row(ct_pB_FE, "sez")$beta,     get_row(ct_pB_FE, "sez")$p_Satt),
    fmt_se( get_row(ct_pB_FE, "sez")$SE),
    n_obs, n_studies, "--", "--", "--", "--"
  ),

  RE = c(
    fmt_est(get_row(ct_pB_RE, "intrcpt")$beta,  get_row(ct_pB_RE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pB_RE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pB_RE, "intrcpt")),
    fmt_est(get_row(ct_pB_RE, "sez")$beta,     get_row(ct_pB_RE, "sez")$p_Satt),
    fmt_se( get_row(ct_pB_RE, "sez")$SE),
    n_obs, n_studies, "--",
    fmt_vc(pB_RE$sigma2), "--",
    fmt_p(lrt_pB_RE$pval)
  ),

  CE = c(
    fmt_est(get_row(ct_pB_CE, "intrcpt")$beta,  get_row(ct_pB_CE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pB_CE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pB_CE, "intrcpt")),
    fmt_est(get_row(ct_pB_CE, "sez")$beta,     get_row(ct_pB_CE, "sez")$p_Satt),
    fmt_se( get_row(ct_pB_CE, "sez")$SE),
    n_obs, n_studies,
    formatC(rho, digits = 3, format = "f"),
    fmt_vc(pB_CE$sigma2), "--", "--"
  ),

  HE = c(
    fmt_est(get_row(ct_pB_HE, "intrcpt")$beta,  get_row(ct_pB_HE, "intrcpt")$p_Satt),
    fmt_se( get_row(ct_pB_HE, "intrcpt")$SE),
    fmt_ci( get_row(ct_pB_HE, "intrcpt")),
    fmt_est(get_row(ct_pB_HE, "sez")$beta,     get_row(ct_pB_HE, "sez")$p_Satt),
    fmt_se( get_row(ct_pB_HE, "sez")$SE