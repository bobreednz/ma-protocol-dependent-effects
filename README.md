# A Protocol for Meta-Analysis with Dependent Effect Sizes

**W. Robert Reed** (Department of Economics and Finance & UCMeta, University of Canterbury; corresponding author, bob.reed@canterbury.ac.nz)  
**Elizabeth Tipton** (School of Education and Social Policy, Northwestern University)  
**František Bartoš** (Faculty of Social and Behavioural Sciences, University of Amsterdam; Institute of Economic Studies, Faculty of Social Sciences, Charles University, Prague)

This repository contains a step-by-step protocol, written in R, for carrying out a meta-analysis when studies report more than one effect size estimate. The protocol is a Quarto document (`Protocol_Combined.qmd`). It explains each analytical decision, sources the R script that carries it out, and displays the resulting tables and figures.

The protocol takes the researcher from the raw data to the final answer. It covers data cleaning and preliminary analysis; the choice of dependence model, comparing five estimators and choosing the correlated and hierarchical effects (CHE) model; three corrections for publication bias (PET-PEESE, step-function selection models, and robust Bayesian meta-analysis, RoBMA); a diagnosis of the selection process that chooses the primary correction; meta-regression and best-practice predictions for all three corrections; and a summary of how sensitive the results are to the choice of correction.

The protocol is illustrated with data from Xue, Reed, and van Aert (2025), a meta-analysis of social capital and economic growth (957 effect sizes from 83 studies). The dataset is only a running example. Every step applies to any meta-analytic database with dependent effect sizes.

**Read the protocol in your browser** at <https://bobreednz.github.io/ma-protocol-dependent-effects/>, or open `Protocol_Combined.html` after downloading.

## What is in this repository

| Type | Files |
|---|---|
| Protocol document | `Protocol_Combined.qmd` (source) and `Protocol_Combined.html` (rendered). The `docs/` folder holds the copy served by GitHub Pages (`protocol.html`) and the landing page (`index.html`). |
| Analysis scripts | 20 `Protocol_*.R` scripts, one per analytical step (listed in the run-order table below), plus `Protocol_Install.R`, which installs the required packages |
| Data | `SCData (20240508).dta` (raw), `SCData_processed.rds` (analysis dataset), `SCData_flagged.dta` (raw data plus quality-check flags) |
| Saved model fits | RoBMA and Bayesian model averaging: `fit_RoBMA_CHE_intercept.rds`, `fit_BMA_CHE_intercept.rds`, `fit_RoBMA_CHE_reg18.rds`, `fit_BMA_CHE_reg18.rds`. Selection models: `sel_3psm_intercept.rds`, `sel_4psm_intercept.rds`, `sel_4psm_reg18_mean.rds`, `sel_4psm_reg18_bp1.rds`, `sel_4psm_reg18_bp2.rds` |
| Outputs | All `Table_*.xlsx` and `Figure_*.png` files, plus `DataCleaning_Report.xlsx` |
| Project files | `MA_Protocol.Rproj`, `.here`, `CITATION.cff` |
| Archive | `archive_HE_RoBMA/` and `archive_superseded/`: scripts and outputs from earlier versions of the protocol, kept for reference and not used by the current document |

The outputs and saved fits are included so that you can render the document and inspect every result immediately. The two RoBMA fits take about 3.5 and 4 hours of MCMC sampling, and the selection-model meta-regression takes 15 to 20 minutes of bootstrapping. Loading the saved files skips all of this.

## Getting started

### 1. Download the files

Click the green **Code** button at the top of this page and choose **Download ZIP** (or clone the repository if you use git). Unzip it to a folder on your computer. Keep all the files together in that one folder.

### 2. Install the software

You need R, RStudio, and Quarto (bundled with recent versions of RStudio). RoBMA also needs JAGS (Just Another Gibbs Sampler), a standalone program that does the MCMC sampling. JAGS is not an R package. Download and install JAGS (version 4.3.1 or newer) from <https://mcmc-jags.sourceforge.io/> before installing RoBMA.

Then install the R packages. The simplest way is to run `Protocol_Install.R`, which installs everything below and prints the versions installed. The same code is:

```r
# Standard packages from CRAN
install.packages(c(
  "here",          # resolves file paths relative to the project root (.Rproj)
  "metafor",       # meta-analytic models
  "clubSandwich",  # cluster-robust (CR2) standard errors
  "tidyverse",     # data manipulation and plotting
  "haven",         # reading Stata .dta files
  "readxl",        # reading Excel output in the display chunks
  "writexl",       # writing Excel output
  "openxlsx",      # writing Excel output with formatting
  "knitr",         # tables and figures in the display chunks
  "future",        # parallel processing (selection-model bootstrap)
  "progressr",     # bootstrap progress display
  "metaselection", # step-function selection models
  "posterior"      # posterior draws for RoBMA predictions
))

# RoBMA 4.1.1 or later, which includes RoBMA.mv(), used to fit the CHE model.
# Install after JAGS, then restart R. (install.packages() also updates an
# older version, such as 4.0.0, if one is already installed.)
install.packages(c("BayesTools", "RoBMA"))
```

### 3. Open the project

Open `MA_Protocol.Rproj` (double-click it, or use File > Open Project in RStudio). This sets the working directory to the project folder, wherever you unzipped it. Every script and the protocol document find their files with `here("filename")`, so there is no path to edit.

### 4. Render the protocol, or reproduce the results

**Render.** Open `Protocol_Combined.qmd` and click **Render**. The document only reads the saved output files, so it renders without running any analysis.

**Reproduce.** Run the scripts in the order below, then render again. Every script after the first reads `SCData_processed.rds`, the dataset produced by script 1, except three that read only the output of earlier scripts: script 11 needs the saved fit from script 10, script 12 needs scripts 5 and 8 to 10, and script 17 needs scripts 12 to 16. The five slowest scripts (9, 10, 15, 16, and 19) load saved results by default, so the whole sequence runs in under an hour on the example data (about half of it for the z-plot, script 11). To re-estimate, set `refit <- TRUE` at the top of the script. The run times in the table are for a full run, with every script refit.

| # | Script | Protocol section | Main output | Full run time |
|:--|:--|:--|:--|:--|
| 1 | `Protocol_DataCleaning.R` | 2 Data cleaning | `SCData_processed.rds`, `DataCleaning_Report.xlsx` | seconds |
| 2 | `Protocol_Variables.R` | 3.1 Variable descriptions | `Table_Variables.xlsx` | seconds |
| 3 | `Protocol_EffectSizes.R` | 3.2 Effect size distribution | `Table_EffectSizes.xlsx`, `Figure_EffectSizes.png` | seconds |
| 4 | `Protocol_StudyWeights.R` | 4.2 Study weights | `Table_StudyWeights.xlsx` | seconds |
| 5 | `Protocol_FiveEstimators.R` | 4 Five models | `Table_FiveEstimators.xlsx` | about 2 minutes |
| 6 | `Protocol_FunnelPlot.R` | 5.1 Funnel plot | `Figure_FunnelPlot.png` | seconds |
| 7 | `Protocol_CHE_PET.R`, `Protocol_CHE_PEESE.R` | 5.2 CHE-PET-PEESE | `Table_CHE_PET.xlsx`, `Table_CHE_PEESE.xlsx` | a few minutes each |
| 8 | `Protocol_CHE_ISCW.R` | 5.2 Robustness check (CHE-ISCW) | `Table_CHE_ISCW.xlsx` | 1 to 3 minutes |
| 9 | `Protocol_SelectionModels.R` | 5.3 Selection models | `Table_SelectionModels.xlsx` | about 3 minutes |
| 10 | `Protocol_RoBMA_CHE_Intercept.R` | 5.4 RoBMA | `Table_RoBMA_CHE_Intercept.xlsx` | about 3.5 hours |
| 11 | `Protocol_Zplot_CHE.R` | 6.1 Z-plot | `Figure_Zplot_CHE.png` | about 5 minutes |
| 12 | `Protocol_Diagnostics.R` | 6.2 Selection weights | `Table_Diagnostics.xlsx` | seconds |
| 13 | `Protocol_CHE_PETPEESE_MetaReg.R` | 7.1 CHE-PET-PEESE meta-regression | `Table_CHE_PETPEESE_MetaReg.xlsx` | 1 to 2 minutes |
| 14 | `Protocol_CHE_PETPEESE_BestPractice.R` | 7.1 Best-practice predictions | `Table_CHE_PETPEESE_BestPractice.xlsx` | under a minute |
| 15 | `Protocol_Selection_MetaReg.R` | 7.2 Selection-model meta-regression | `Table_Selection_MetaReg.xlsx` | 15 to 20 minutes |
| 16 | `Protocol_RoBMA_CHE_MetaReg.R` | 7.3 RoBMA meta-regression | `Table_RoBMA_CHE_MetaReg.xlsx` | about 4 hours |
| 17 | `Protocol_Sensitivity.R` | 8 Sensitivity summary | `Table_Sensitivity.xlsx` | seconds |
| 18 | `Protocol_RhoSensitivity.R` | Appendix A: choice of rho | `Table_RhoSensitivity.xlsx`, two figures | about 6 minutes |
| 19 | `Protocol_Outliers.R` | Appendix B: outliers | `Table_Jackknife_CHE.xlsx`, two figures | about 1.5 hours |

Run times are for the example data on a standard desktop computer; with other data they vary with the numbers of estimates, studies, and moderators. The selection-model bootstraps run in parallel and take longer with fewer processor cores.

## Software versions

The results in the protocol were produced with RoBMA 4.1.1 and BayesTools 0.4.0 from CRAN, and with metaselection 0.1.5. Later versions of these packages may give slightly different results. Run `packageVersion("RoBMA")` to see which version you have installed.

## Using this protocol with your own data

The protocol is written as a template. Replace the data file, adjust the variable names in the scripts (the key inputs are an effect size, its standard error, and a study identifier), and work through the same sequence of decisions. Start with `Protocol_DataCleaning.R`, which writes the analysis dataset that every later script reads. Then set `refit <- TRUE` at the top of scripts 9, 10, 15, 16, and 19. These scripts load saved results by default, so without this change they would report the results for the example data.

Each section of the protocol has the same structure. "What to look for" explains how to read the output for any dataset. "Decision point" sets out a choice the researcher must make before going on. "In this example" reports what the output shows for the social capital data. The three corrections for publication bias are always estimated in the same order, so the document works the same way whichever correction your own diagnostics favor. Section 5.4 explains how to use RoBMA with effect sizes other than Fisher's z.

## Citation

If you use this protocol, please cite the repository (GitHub's **Cite this repository** button uses `CITATION.cff`). For the illustrative data, cite:

Xue, X., W. R. Reed, and R. C. M. van Aert. 2025. "Social Capital and Economic Growth: A Meta-Analysis." *Journal of Economic Surveys* 39(4): 1395-1432.

## License

The protocol, scripts, and outputs are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0); see `LICENSE`. You may reuse and adapt them, provided you give appropriate credit. The illustrative dataset comes from Xue, Reed, and van Aert (2025); please cite that article when you use it.

## Questions

Contact W. Robert Reed at bob.reed@canterbury.ac.nz.
