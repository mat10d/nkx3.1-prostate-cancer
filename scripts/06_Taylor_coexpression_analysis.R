################################################################################
# Taylor 2010 - NKX3.1 and HSPA9 Co-expression Analysis
################################################################################
#
# Project: NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer
#          Initiation
# PI: Cory Abate-Shen
# Publication: Cancer Discovery (2021)
# DOI: https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307
#
# Description:
# This script analyzes the correlation between NKX3.1 and HSPA9B expression
# in the Taylor et al. 2010 prostate cancer dataset (GSE21035). The dataset
# contains microarray data, and HSPA9B has three probe sets that are analyzed
# separately.
#
# Reference:
# Taylor BS et al. (2010) Integrative genomic profiling of human prostate cancer.
# Cancer Cell 18(1):11-22.
#
################################################################################

# Load required libraries
library(Biobase)
library(dplyr)
library(ggpubr)

################################################################################
# SECTION 1: Data Loading
################################################################################

# Load Taylor 2010 processed data (GSE21035)
taylor <- readRDS("../data/taylor-2010-processed-data.rds")

# Extract expression matrix
exprs_taylor <- taylor$exprs

# Extract sample metadata
columndata_taylor <- taylor$metadata

# Verify sample matching
stopifnot(identical(colnames(exprs_taylor), rownames(columndata_taylor)))

# Create output directory
dir.create("../results/Co-expression", recursive = TRUE, showWarnings = FALSE)

################################################################################
# SECTION 2: Select Gene Probes
################################################################################

# Select probes for NKX3-1 and HSPA9B
# Note: HSPA9B has three different probe sets on the array
select_genes <- taylor$probes[taylor$probes$GENE_SYMBOL %in% c("NKX3-1", "HSPA9B"), ]

# Extract expression data for selected genes
exprs_taylor_reduced <- as.data.frame(t(exprs_taylor[select_genes$ID, ]))

# Assign informative column names
# Three HSPA9B probes: A_14_P128569, A_14_P201971, A_16_P37383218
colnames(exprs_taylor_reduced) <- c("NKX3.1", "HSPA9Ba", "HSPA9Bb", "HSPA9Bc")

################################################################################
# SECTION 3: Correlation Analysis - HSPA9B Probe A (A_14_P128569)
################################################################################

# Perform Pearson correlation test
res_a <- cor.test(exprs_taylor_reduced$NKX3.1, exprs_taylor_reduced$HSPA9Ba,
                  method = "pearson")

# Save correlation results
sink('../results/Co-expression/Taylor_correlation_A_14_P128569_statistics.txt')
cat("=== Taylor 2010: NKX3.1 and HSPA9B Correlation ===\n")
cat("Probe: A_14_P128569\n\n")
cat("Number of samples:", nrow(exprs_taylor_reduced), "\n\n")
cat("Pearson Correlation Test\n")
cat("------------------------\n\n")
print(res_a)
sink()

# Generate scatter plot
p_a <- ggscatter(
  exprs_taylor_reduced,
  x = "NKX3.1",
  y = "HSPA9Ba",
  fullrange = TRUE,
  add = "reg.line",
  conf.int = FALSE,
  size = 0.5,
  cor.coef.coord = c(-0.28, 0.375),
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_bw()
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black")
  ) +
  ggplot2::annotate("text", x = -0.18, y = 0.35,
                    label = paste("n =", nrow(exprs_taylor_reduced)), size = 4) +
  xlab("NKX3.1 Expression") +
  ylab("HSPA9B Expression (Probe A_14_P128569)") +
  ggtitle("Taylor 2010: NKX3.1 vs HSPA9B (Probe A)")

ggsave("../results/Co-expression/Taylor_correlation_A_14_P128569.pdf",
       plot = p_a, height = 5, width = 5)

################################################################################
# SECTION 4: Correlation Analysis - HSPA9B Probe B (A_14_P201971)
################################################################################

# Perform Pearson correlation test
res_b <- cor.test(exprs_taylor_reduced$NKX3.1, exprs_taylor_reduced$HSPA9Bb,
                  method = "pearson")

# Save correlation results
sink('../results/Co-expression/Taylor_correlation_A_14_P201971_statistics.txt')
cat("=== Taylor 2010: NKX3.1 and HSPA9B Correlation ===\n")
cat("Probe: A_14_P201971\n\n")
cat("Number of samples:", nrow(exprs_taylor_reduced), "\n\n")
cat("Pearson Correlation Test\n")
cat("------------------------\n\n")
print(res_b)
sink()

# Generate scatter plot
p_b <- ggscatter(
  exprs_taylor_reduced,
  x = "NKX3.1",
  y = "HSPA9Bb",
  fullrange = TRUE,
  add = "reg.line",
  conf.int = FALSE,
  size = 0.5,
  cor.coef.coord = c(-0.28, 0.375),
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_bw()
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black")
  ) +
  ggplot2::annotate("text", x = -0.18, y = 0.35,
                    label = paste("n =", nrow(exprs_taylor_reduced)), size = 4) +
  xlab("NKX3.1 Expression") +
  ylab("HSPA9B Expression (Probe A_14_P201971)") +
  ggtitle("Taylor 2010: NKX3.1 vs HSPA9B (Probe B)")

ggsave("../results/Co-expression/Taylor_correlation_A_14_P201971.pdf",
       plot = p_b, height = 5, width = 5)

################################################################################
# SECTION 5: Correlation Analysis - HSPA9B Probe C (A_16_P37383218)
################################################################################

# Perform Pearson correlation test
res_c <- cor.test(exprs_taylor_reduced$NKX3.1, exprs_taylor_reduced$HSPA9Bc,
                  method = "pearson")

# Save correlation results
sink('../results/Co-expression/Taylor_correlation_A_16_P37383218_statistics.txt')
cat("=== Taylor 2010: NKX3.1 and HSPA9B Correlation ===\n")
cat("Probe: A_16_P37383218\n\n")
cat("Number of samples:", nrow(exprs_taylor_reduced), "\n\n")
cat("Pearson Correlation Test\n")
cat("------------------------\n\n")
print(res_c)
sink()

# Generate scatter plot
p_c <- ggscatter(
  exprs_taylor_reduced,
  x = "NKX3.1",
  y = "HSPA9Bc",
  fullrange = TRUE,
  add = "reg.line",
  conf.int = FALSE,
  size = 0.5,
  cor.coef.coord = c(-0.28, 0.375),
  cor.coef = TRUE,
  cor.method = "pearson",
  ggtheme = theme_bw()
) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black")
  ) +
  ggplot2::annotate("text", x = -0.18, y = 0.35,
                    label = paste("n =", nrow(exprs_taylor_reduced)), size = 4) +
  xlab("NKX3.1 Expression") +
  ylab("HSPA9B Expression (Probe A_16_P37383218)") +
  ggtitle("Taylor 2010: NKX3.1 vs HSPA9B (Probe C)")

ggsave("../results/Co-expression/Taylor_correlation_A_16_P37383218.pdf",
       plot = p_c, height = 5, width = 5)

################################################################################
# SECTION 6: Save Summary Statistics
################################################################################

# Create summary data frame with correlation statistics for all probes
summary_df <- data.frame(
  Dataset = c("Taylor 2010 (Probe A)", "Taylor 2010 (Probe B)", "Taylor 2010 (Probe C)"),
  Probe_ID = c("A_14_P128569", "A_14_P201971", "A_16_P37383218"),
  N_Samples = rep(nrow(exprs_taylor_reduced), 3),
  Pearson_R = c(res_a$estimate, res_b$estimate, res_c$estimate),
  P_value = c(res_a$p.value, res_b$p.value, res_c$p.value),
  CI_Lower = c(res_a$conf.int[1], res_b$conf.int[1], res_c$conf.int[1]),
  CI_Upper = c(res_a$conf.int[2], res_b$conf.int[2], res_c$conf.int[2])
)

write.csv(summary_df, "../results/Co-expression/Taylor_correlation_summary.csv",
          row.names = FALSE)

################################################################################
# Analysis Complete
################################################################################

cat("\nTaylor 2010 co-expression analysis completed successfully.\n")
cat("Results saved to: ../results/Co-expression/\n\n")
cat("Correlation Summary:\n")
cat("  Probe A (A_14_P128569):   R =", round(res_a$estimate, 3),
    ", p =", format(res_a$p.value, scientific = TRUE), "\n")
cat("  Probe B (A_14_P201971):   R =", round(res_b$estimate, 3),
    ", p =", format(res_b$p.value, scientific = TRUE), "\n")
cat("  Probe C (A_16_P37383218): R =", round(res_c$estimate, 3),
    ", p =", format(res_c$p.value, scientific = TRUE), "\n")
