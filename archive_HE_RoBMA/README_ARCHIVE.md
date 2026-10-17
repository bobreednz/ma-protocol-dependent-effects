# Archived: HE RoBMA analysis (retired 2026-10-14)

These scripts and outputs belong to the pre-October 2026 version of the protocol, which fitted multilevel RoBMA with an HE (hierarchical effects) dependence structure (RoBMA 4.0.0, `cluster = newid`), a BIC moderator pre-screen (Table 9), and an Appendix C robustness check with an alternative moderator set.

They were retired when RoBMA 4.1 made CHE models available (`RoBMA.mv()`, `BMA.mv()`), and the protocol was restructured: three publication-bias corrections at the intercept-only stage, diagnosis and choice of a primary method, then meta-regression and best-practice predictions for all three methods.

Kept for reference only. They are no longer sourced by `Protocol_Combined.qmd`.

Note: several outputs here were made with sampling variances 1/(df - 1); the current protocol uses sez^2 throughout.

| Script | Outputs |
|---|---|
| Protocol_Figure3.R | Figure3_Zplot.png, RoBMA_intercept_only.rds |
| Protocol_Table8.R | Table8_RoBMA_Comparison.xlsx |
| Protocol_Table9.R | Table9_BIC_ModeratorSelection.xlsx |
| Protocol_Table10.R | Table10_RoBMA_MetaRegression.xlsx, RoBMA_metaregression.rds |
| Protocol_Table11.R | Table11_RoBMA_BestPractice.xlsx |
| Protocol_Table10_FYI.R | Table10_FYI_RoBMA.xlsx, RoBMA_metaregression_FYI.rds |
| Protocol_Table11_FYI.R | Table11_FYI_BestPractice.xlsx |
| Protocol_Table12.R | Table12_ComparisonSummary.xlsx |
