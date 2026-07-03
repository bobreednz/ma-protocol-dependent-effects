# A Protocol for Meta-Analysis with Dependent Effect Sizes

**W. Robert Reed** (Department of Economics and Finance & UCMeta, University of Canterbury) — bob.reed@canterbury.ac.nz

This repository contains a complete, step-by-step protocol for conducting a meta-analysis when effect size estimates are statistically dependent (multiple estimates per study). The protocol is written as a Quarto document (`Protocol_Combined.qmd`) that walks through every analytical decision, sources the R script for each step, and displays the resulting tables and figures. It covers data cleaning, preliminary analysis, the choice of the within-study correlation (rho), five meta-analytic estimators (FE, RE, CE, HE, and the preferred CHE model), publication bias analysis using PET-PEESE plus robustness checks (CHE-ISCW weighting and step-function selection models), meta-regression combining Bayesian Model Averaging with CHE, best-practice effect size predictions, and a full Robust Bayesian Meta-Analysis (RoBMA) extension with a concluding comparison of the PET-PEESE and RoBMA corrections.

The protocol is illustrated with data from Xue, Reed, and van Aert (2024), a meta-analysis of the effect of social capital on economic growth (957 effect sizes from 83 studies), but every step applies to any meta-analytic database with dependent effect sizes.

**Read the rendered protocol in your browser** at the GitHub Pages site for this repository, or open `Protocol_Combined.html` after downloading.

## What is in this repository

| Type | Files |
|---|---|
| Protocol document | `Protocol_Combined.qmd` (source), `Protocol_Combined.html` (rendered) |
| Analysis scripts | 18 `Protocol_*.R` scripts, one per analytical step |
| Data | `SCData (20240508).dta` (raw), `SCData_processed.rds` (analysis dataset), `SCData_flagged.dta` (raw data plus quality-check flags) |
| Pre-computed outputs | All `Table*.xlsx` and `Figure*.png` files, plus `DataCleaning_Report.xlsx` |
| Fitted RoBMA models | `RoBMA_intercept_only.rds`, `RoBMA_metaregression.rds` |

The pre-computed outputs and fitted RoBMA models are included so that you can render the full document and inspect every result **immediately**, without first running the computationally expensive steps. The two RoBMA model fits alone take roughly 5 and 4 hours of MCMC sampling; loading the saved `.rds` files skips that entirely.

## Getting started

### 1. Download the files

Click the green **Code** button at the top of this page and choose **Download ZIP** (or clone the repository if you use git). Unzip to a folder on your computer, for example `C:/MA_PROTOCOL`.

### 2. Install the software

You need R (version 4.2 or later recommended), RStudio, and Quarto (bundled with recent RStudio versions). Then install the required packages in R:

```r
install.packages(c(
  "metafor",       # meta-analytic models
  "clubSandwich",  # cluster-robust (CR2) standard errors
  "tidyverse",     # data manipulation and plotting
  "haven",         # reading Stata .dta files
  "writexl",       # writing Excel output
  "openxlsx",      # writing Excel output with formatting
  "kableExtra",    # table formatting in display chunks
  "BMS",           # Bayesian Model Averaging
  "RoBMA",         # Robust Bayesian meta-analysis
  "future",        # parallel processing
  "progressr"      # bootstrap progress display
))

# metaselection is on GitHub only (not CRAN):
# install.packages("remotes")
remotes::install_github("jepusto/metaselection")
```

Note: `RoBMA` requires the JAGS sampler. If `install.packages("RoBMA")` reports a JAGS problem, install JAGS first from https://mcmc-jags.sourceforge.io and then reinstall RoBMA.

### 3. Point the scripts at your folder

Every script begins with a `setwd()` line that sets the working directory:

```r
setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")
```

Replace that path with the folder where you unzipped the files. The fastest way is RStudio's **Edit > Find in Files** (Ctrl+Shift+F): search for `C:/PROTOCOL OF MAs WITH DEPENDENT DATA`, set the directory to your unzipped folder, and replace all occurrences with your own path. The same path appears once in the `setup` chunk of `Protocol_Combined.qmd`.

### 4. Open the protocol and work through it

Open `Protocol_Combined.qmd` in RStudio. You can use it two ways.

**Read and render.** Click **Render**. Because all output files are included in the repository, the document renders in a minute or two without running any analysis, reproducing the full HTML with every table and figure.

**Reproduce the results.** Work through the document section by section. Each section has a "Run the script" chunk (`source("Protocol_XXX.R")`). Run these in order to regenerate every output file from the raw data, then re-render the document. Approximate run times on a standard desktop:

| Script | Step | Run time |
|---|---|---|
| `Protocol_DataCleaning.R` | Data cleaning and quality checks | seconds |
| `Protocol_Table1.R` | Variable descriptions | seconds |
| `Protocol_Table2_Figure1.R` | Effect size descriptives and histogram | seconds |
| `Protocol_Table3.R` | Study weights (FE, RE, HE) | ~1 min |
| `Protocol_RhoSensitivity.R` | Profile likelihood over rho | ~6 min |
| `Protocol_Table4.R` | The five estimators | ~2 min |
| `Protocol_Figure2.R` | Funnel plot | seconds |
| `Protocol_Table5.R` | FAT-PET regressions | ~4 min |
| `Protocol_Table5S.R` | PEESE sensitivity | ~4 min |
| `Protocol_Table6.R` | BMA + CHE meta-regression | ~5 min |
| `Protocol_Table7.R` | Best-practice predictions (CHE) | ~40 sec |
| `Protocol_AppendixPubBias.R` | Robustness checks (ISCW, 3PSM/4PSM) | ~5-30 min (bootstrap; uses parallel cores) |
| `Protocol_Figure3.R` | Intercept-only RoBMA + z-plot | **~5 hours** |
| `Protocol_Table8.R` | Comparison of corrections | seconds |
| `Protocol_Table9.R` | BIC moderator pre-selection | ~15 min |
| `Protocol_Table10.R` | RoBMA meta-regression | **~4 hours** |
| `Protocol_Table11.R` | RoBMA best-practice predictions | ~1 min |
| `Protocol_Outliers.R` | Influence diagnostics + jackknife | **~2 hours** |

The three long-running scripts save their fitted models (`.rds` files), which downstream scripts load instead of refitting. Since those `.rds` files are included in the repository, you can reproduce Tables 8 and 11 without running the 5-hour and 4-hour fits — only run `Protocol_Figure3.R` and `Protocol_Table10.R` if you want to verify the MCMC estimation itself.

## Using this protocol with your own data

The protocol is written to be adapted. Replace the data file, adjust the variable names in the scripts (the key inputs are an effect size, its standard error, and a study identifier), and work through the same sequence of decisions. Section 1 of the protocol describes what changes when your effect size metric differs from the partial correlations used here.

## Citation

If you use this protocol, please cite the repository, and for the illustrative data:

Xue, H., W. Reed, and R. C. J. van Aert. 2024. "The Effects of Social Capital on Economic Growth: A Meta-Analysis." *Journal of Economic Surveys* 38: 1234-1268.

## Questions

Contact W. Robert Reed at bob.reed@canterbury.ac.nz.
