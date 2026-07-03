setwd("C:/PROTOCOL OF MAs WITH DEPENDENT DATA")

# Record start time -- printed at the end so run time can be documented.
start_time <- proc.time()

# ============================================================
# Install and load required packages
# ============================================================

packages <- c("tidyverse", "writexl")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)  # ggplot2, dplyr, etc.
library(writexl)    # write_xlsx() for Excel export

# ============================================================
# Load data
# ============================================================

DT <- readRDS("SCData_processed.rds")

# ============================================================
# TABLE 2: Descriptive statistics for effect sizes
#
# We report Fisher's z (the primary effect size measure) and
# PCC (partial correlation coefficient) side by side. PCC is
# included because it is more interpretable to substantive
# readers: a PCC of 0.07 corresponds to a "small" effect by
# the Doucouliagos (2011) thresholds. Fisher's z and PCC are
# nearly identical for small values but diverge for larger ones.
#
# We report mean, median, standard deviation, minimum,
# maximum, and observation count. This is intentionally more
# concise than the paper's Table 3, which includes t-statistics,
# degrees of freedom, and a long list of percentiles.
# ============================================================

# Helper: compute the five summary statistics for a vector
summarise_var <- function(x) {
  c(mean(x), median(x), sd(x), min(x), max(x))
}

# Format to 3 decimal places
fmt3 <- function(x) formatC(x, digits = 3, format = "f")

table2 <- tibble(
  Statistic      = c("Mean", "Median", "Std. Dev.", "Minimum",
                     "Maximum", "Observations"),
  `t-statistics` = c(sapply(summarise_var(DT$tstat), fmt3),
                     as.character(nrow(DT))),
  `df`           = c(sapply(summarise_var(DT$df),    fmt3),
                     as.character(nrow(DT))),
  `PCC values`   = c(sapply(summarise_var(DT$pcc),   fmt3),
                     as.character(nrow(DT))),
  `z values`     = c(sapply(summarise_var(DT$z),     fmt3),
                     as.character(nrow(DT)))
)

write_xlsx(table2, "Table2_Protocol.xlsx")
cat("Table 2 saved to Table2_Protocol.xlsx\n")

# ============================================================
# FIGURE 1: Distribution of Fisher's z values
#
# The histogram shows the distribution of effect sizes across
# all 957 estimates. Two sets of reference lines are included:
#
#   Dashed lines: Fisher's z = +/-0.07, the approximate
#     equivalent of Doucouliagos' (2011) "small" PCC threshold
#     (Table 3, column 1, "All estimates").
#     (For small values, Fisher's z and PCC are nearly equal,
#     so this correspondence is close but not exact.)
#
#   Solid line: the mean Fisher's z (computed dynamically from
#     the data), giving the reader an immediate visual anchor
#     for the central tendency of the distribution.
#
# The inset table reports the share of estimates whose
# t-statistic falls below -2, between -2 and 2, or above 2.
# This approximates two-tailed significance at the 5 percent
# level. The pronounced asymmetry between the positive and
# negative tails (45.9% vs. 5.2%) is a first visual signal of
# publication bias and motivates the funnel plots and FAT-PET
# tests that follow.
# ============================================================

# --- Significance breakdown (based on t-statistics) ---
pct_neg <- round(mean(DT$tstat < -2.00) * 100, 1)
pct_mid <- round(mean(DT$tstat >= -2.00 & DT$tstat <= 2.00) * 100, 1)
pct_pos <- round(mean(DT$tstat >  2.00) * 100, 1)

mean_z <- mean(DT$z)

# --- Build the histogram ---
# after_stat(count / sum(count) * 100) converts counts to
# percentages on the y-axis.
# binwidth = 0.05 gives a smooth but detailed picture. Adjust
# if the histogram looks too jagged (increase) or too coarse
# (decrease) for a different dataset.
# boundary = 0 aligns bins to round numbers on the x-axis.

fig1 <- ggplot(DT, aes(x = z)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count) * 100)),
    binwidth = 0.05,
    fill     = "gray60",
    color    = "white",
    boundary = 0
  ) +
  # Dashed lines at +/-0.07 (Doucouliagos 2011 "small" threshold)
  geom_vline(xintercept = c(-0.07, 0.07),
             linetype = "dashed", linewidth = 0.6) +
  # Solid line at the mean Fisher's z
  geom_vline(xintercept = mean_z,
             linetype = "solid", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(-1, 2, by = 0.5)) +
  labs(x = "Fisher's z", y = "Percent") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# --- Add the inset significance table ---
# The inset is drawn using annotate() calls: one rectangle for
# the box background, then separate text and segment calls for
# the header, divider line, and data rows.
#
# The coordinates below are calibrated for this dataset, where
# the y-axis reaches roughly 10 percent and the x-axis spans
# approximately -1 to 2. If you apply this script to a dataset
# with a very different scale, adjust xmin/xmax/ymin/ymax.

fig1 <- fig1 +
  # White background box
  annotate("rect",
           xmin = 0.80, xmax = 1.97,
           ymin = 5.0,  ymax = 10.3,
           fill = "white", color = "black", linewidth = 0.4) +
  # Title row: spans the full width of the box, centred
  annotate("text",
           x = (0.80 + 1.97) / 2, y = 9.8,
           hjust = 0.5, size = 2.8, fontface = "bold",
           label = "Distribution of t-statistics") +
  # Horizontal line below title
  annotate("segment",
           x = 0.82, xend = 1.94,
           y = 9.2,  yend = 9.2,
           linewidth = 0.3) +
  # Row 1: negative significant
  annotate("text", x = 0.87, y = 8.5, hjust = 0, size = 2.8,
           label = "t < -2.00") +
  annotate("text", x = 1.91, y = 8.5, hjust = 1, size = 2.8,
           label = paste0(pct_neg, "%")) +
  # Row 2: not significant
  annotate("text", x = 0.87, y = 7.5, hjust = 0, size = 2.8,
           label = "-2.00 <= t <= 2.00") +
  annotate("text", x = 1.91, y = 7.5, hjust = 1, size = 2.8,
           label = paste0(pct_mid, "%")) +
  # Row 3: positive significant
  annotate("text", x = 0.87, y = 6.5, hjust = 0, size = 2.8,
           label = "t > 2.00") +
  annotate("text", x = 1.91, y = 6.5, hjust = 1, size = 2.8,
           label = paste0(pct_pos, "%"))

# --- Save the figure ---
# width = 6, height = 4 inches at 300 dpi gives a figure
# suitable for most journal and document formats.

ggsave("Figure1_Protocol.png", fig1,
       width = 6, height = 4, dpi = 300)

cat("Figure 1 saved to Figure1_Protocol.png\n")

# ── Run time ──────────────────────────────────────────────────────────────────────────────
elapsed <- pro