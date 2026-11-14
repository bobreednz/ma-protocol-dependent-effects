# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("here", "tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(here)       # resolves file paths relative to the project root (.Rproj)
library(tidyverse)
library(writexl)

# ============================================================
# Load data
# ============================================================

DT <- readRDS(here("SCData_processed.rds"))

n_studies <- length(unique(DT$newid))

# ============================================================
# Helper: count distinct studies where a binary variable = 1
# ============================================================

# For binary variables, "Studies" reports how many unique
# studies (newid) have at least one effect size where the
# variable equals 1. This allows a study to appear under
# more than one category (e.g., a study contributing both
# cognitive and structural SC estimates is counted in both).

n_studies_eq1 <- function(var) {
  length(unique(DT$newid[DT[[var]] == 1]))
}

# ============================================================
# Build the table
#
# Structure mirrors Table 2 in Xue, Reed, and van Aert (2025)
# with one addition: a "Studies" column inserted before "Mean".
# For binary variables, Studies = number of distinct studies
# with at least one observation equal to 1.
# For continuous variables (PubYear, NumberSCVars), Studies
# is left blank because the concept does not apply.
#
# Reference categories (omitted in the regression) are marked
# with an asterisk in the Variable column, matching the paper.
# ============================================================

rows <- tribble(
  ~Group, ~Variable, ~Description, ~Studies, ~Mean, ~Min, ~Max,

  # ---- Economic Growth Variable ----
  "ECONOMIC GROWTH VARIABLE",
  "DV_GrowthRate",
  "=1, if dep. variable is GDP growth rate",
  n_studies_eq1("DV_GrowthRate"),
  mean(DT$DV_GrowthRate), min(DT$DV_GrowthRate), max(DT$DV_GrowthRate),

  "ECONOMIC GROWTH VARIABLE",
  "DV_GDPLevel*",
  "=1, if dep. variable is level of income (GDP)",
  n_studies_eq1("DV_GDPLevel"),
  mean(DT$DV_GDPLevel), min(DT$DV_GDPLevel), max(DT$DV_GDPLevel),

  # ---- Social Capital Variables ----
  "SOCIAL CAPITAL VARIABLES",
  "SC1_Cognitive",
  "=1, if social capital is cognitive",
  n_studies_eq1("SC1_Cognitive"),
  mean(DT$SC1_Cognitive), min(DT$SC1_Cognitive), max(DT$SC1_Cognitive),

  "SOCIAL CAPITAL VARIABLES",
  "SC1_Structural",
  "=1, if social capital is structural",
  n_studies_eq1("SC1_Structural"),
  mean(DT$SC1_Structural), min(DT$SC1_Structural), max(DT$SC1_Structural),

  "SOCIAL CAPITAL VARIABLES",
  "SC1_Other*",
  "=1, if social capital is other type",
  n_studies_eq1("SC1_Other"),
  mean(DT$SC1_Other), min(DT$SC1_Other), max(DT$SC1_Other),

  # ---- Study Characteristics ----
  "STUDY CHARACTERISTICS",
  "PubYear",
  "Year study was published",
  NA_integer_,
  mean(DT$PubYear), min(DT$PubYear), max(DT$PubYear),

  "STUDY CHARACTERISTICS",
  "Published",
  "=1, if study is a published journal article",
  n_studies_eq1("Published"),
  mean(DT$Published), min(DT$Published), max(DT$Published),

  # ---- Specification ----
  "SPECIFICATION",
  "LaggedDV",
  "=1, if lagged DV included in equation",
  n_studies_eq1("LaggedDV"),
  mean(DT$LaggedDV), min(DT$LaggedDV), max(DT$LaggedDV),

  "SPECIFICATION",
  "LaggedSC",
  "=1, if lagged SC variable(s) included in equation",
  n_studies_eq1("LaggedSC"),
  mean(DT$LaggedSC), min(DT$LaggedSC), max(DT$LaggedSC),

  "SPECIFICATION",
  "NumberSCVars",
  "Number of SC variables included in equation",
  NA_integer_,
  mean(DT$NumberSCVars), min(DT$NumberSCVars), max(DT$NumberSCVars),

  # ---- Endogeneity ----
  "ENDOGENEITY",
  "Endog_IV",
  "=1, if instrumental variable estimator used",
  n_studies_eq1("Endog_IV"),
  mean(DT$Endog_IV), min(DT$Endog_IV), max(DT$Endog_IV),

  "ENDOGENEITY",
  "Endog_FE",
  "=1, if fixed effects included in equation",
  n_studies_eq1("Endog_FE"),
  mean(DT$Endog_FE), min(DT$Endog_FE), max(DT$Endog_FE),

  "ENDOGENEITY",
  "NoEndogeneity*",
  "=1, if estimation did not address endogeneity",
  n_studies_eq1("NoEndogeneity"),
  mean(DT$NoEndogeneity), min(DT$NoEndogeneity), max(DT$NoEndogeneity),

  "ENDOGENEITY",
  "PanelData",
  "=1, if data are panel data",
  n_studies_eq1("PanelData"),
  mean(DT$PanelData), min(DT$PanelData), max(DT$PanelData),

  # ---- Spatial Data Characteristics ----
  "SPATIAL DATA CHARACTERISTICS",
  "CityLevel",
  "=1, if data are city level",
  n_studies_eq1("CityLevel"),
  mean(DT$CityLevel), min(DT$CityLevel), max(DT$CityLevel),

  "SPATIAL DATA CHARACTERISTICS",
  "RegionLevel",
  "=1, if data are regional level",
  n_studies_eq1("RegionLevel"),
  mean(DT$RegionLevel), min(DT$RegionLevel), max(DT$RegionLevel),

  "SPATIAL DATA CHARACTERISTICS",
  "CountryLevel",
  "=1, if data are country level",
  n_studies_eq1("CountryLevel"),
  mean(DT$CountryLevel), min(DT$CountryLevel), max(DT$CountryLevel),

  "SPATIAL DATA CHARACTERISTICS",
  "OtherLevel*",
  "=1, if data are other level",
  n_studies_eq1("OtherLevel"),
  mean(DT$OtherLevel), min(DT$OtherLevel), max(DT$OtherLevel),

  # ---- Countries ----
  "COUNTRIES",
  "Reg_OECDEurope",
  "=1, if countries are in OECD or Europe",
  n_studies_eq1("Reg_OECDEurope"),
  mean(DT$Reg_OECDEurope), min(DT$Reg_OECDEurope), max(DT$Reg_OECDEurope),

  "COUNTRIES",
  "Reg_US",
  "=1, if countries are in US",
  n_studies_eq1("Reg_US"),
  mean(DT$Reg_US), min(DT$Reg_US), max(DT$Reg_US),

  "COUNTRIES",
  "Reg_Africa",
  "=1, if countries are in Africa",
  n_studies_eq1("Reg_Africa"),
  mean(DT$Reg_Africa), min(DT$Reg_Africa), max(DT$Reg_Africa),

  "COUNTRIES",
  "Reg_Asia",
  "=1, if countries are in Asia",
  n_studies_eq1("Reg_Asia"),
  mean(DT$Reg_Asia), min(DT$Reg_Asia), max(DT$Reg_Asia),

  "COUNTRIES",
  "Reg_Other*",
  "=1, if countries outside the above regions",
  n_studies_eq1("Reg_Other"),
  mean(DT$Reg_Other), 0, 1
)

# ============================================================
# Format numeric columns for display
# ============================================================

# Mean: 3 decimal places for binary vars; 1 decimal for
# continuous vars (PubYear, NumberSCVars)
# Min/Max: computed from the data for every variable, so a binary
# variable with a stray value (e.g. 2) shows up in the table.
# Integers for binary vars; no decimals for others

table_out <- rows %>%
  mutate(
    Studies = if_else(is.na(Studies), "", as.character(Studies)),
    Mean    = case_when(
      Variable %in% c("PubYear", "NumberSCVars") ~
        formatC(Mean, digits = 1, format = "f"),
      TRUE ~
        formatC(Mean, digits = 3, format = "f")
    ),
    Min = case_when(
      Variable == "PubYear"       ~ as.character(as.integer(Min)),
      Variable == "NumberSCVars"  ~ as.character(as.integer(Min)),
      TRUE                        ~ as.character(as.integer(Min))
    ),
    Max = case_when(
      Variable == "PubYear"       ~ as.character(as.integer(Max)),
      Variable == "NumberSCVars"  ~ as.character(as.integer(Max)),
      TRUE                        ~ as.character(as.integer(Max))
    )
  )

# ============================================================
# Add note row explaining the asterisk convention
# ============================================================

# A blank row followed by the note keeps it visually separate
# from the data rows. The note text goes in the Group column;
# all other columns are left empty.

note_rows <- tibble(
  Group    = c("",
               "Note: Asterisks (*) indicate benchmark (reference) categories omitted in the regression analysis."),
  Variable = c("", ""),
  Description = c("", ""),
  Studies  = c("", ""),
  Mean     = c("", ""),
  Min      = c("", ""),
  Max      = c("", "")
)

table_out <- bind_rows(table_out, note_rows)

# ============================================================
# Export to Excel
# ============================================================

write_xlsx(table_out, here("Table_Variables.xlsx"))

cat("Done. Output saved to Table_Variables.xlsx\n")
cat("Total studies:", n_studies, "\n")

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- proc.time() - start_time
cat(sprintf("\nTotal run time: %.1f seconds (%.1f minutes).\n",
            elapsed["elapsed"], elapsed["elapsed"] / 60))
