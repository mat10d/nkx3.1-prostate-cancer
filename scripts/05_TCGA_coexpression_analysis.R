################################################################################
# TCGA PRAD - NKX3.1 and HSPA9 Co-expression Analysis
################################################################################
#
# Project: NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer
#          Initiation
# PI: Cory Abate-Shen
# Publication: Cancer Discovery (2021)
# DOI: https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307
#
# Description:
# This script analyzes the correlation between NKX3.1 and HSPA9 expression
# in TCGA prostate adenocarcinoma (PRAD) samples. It generates scatter plots
# showing the relationship between these genes using both raw counts and
# variance-stabilized normalized expression values.
#
################################################################################

# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(data.table)

################################################################################
# SECTION 1: Data Loading
################################################################################

# Load TCGA PRAD gene expression data
load("../data/TCGA-PRAD-GeneExpression-Counts.rda")

prad_data <- data

# Extract raw counts matrix
raw_counts_tcga <- assay(prad_data)

# Extract column metadata
TCGA_PRAD_DATA <- as.data.frame(colData(prad_data))

# Verify sample matching
stopifnot(identical(colnames(raw_counts_tcga), rownames(TCGA_PRAD_DATA)))

# Create output directory
dir.create("../results/Co-expression", recursive = TRUE, showWarnings = FALSE)

################################################################################
# SECTION 2: Variance Stabilizing Transformation
################################################################################

# Apply variance stabilizing transformation for visualization
# This transformation normalizes count data and stabilizes variance across
# the mean, making it suitable for correlation analysis
vsd <- vst(raw_counts_tcga, blind = FALSE)

################################################################################
# SECTION 3: Extract NKX3.1 and HSPA9 Expression
################################################################################

# Gene IDs:
# ENSG00000167034 = HSPA9 (Heat Shock Protein Family A Member 9)
# ENSG00000113013 = NKX3.1 (NK3 Homeobox 1)

# Extract raw counts for HSPA9 and NKX3.1
TCGA_reduced <- as.data.frame(t(raw_counts_tcga[
  rownames(raw_counts_tcga) %in% c("ENSG00000167034", "ENSG00000113013"), ]))

# Extract variance-stabilized expression for HSPA9 and NKX3.1
vsd_reduced <- as.data.frame(t(vsd[
  rownames(vsd) %in% c("ENSG00000167034", "ENSG00000113013"), ]))

# Rename columns to gene symbols
colnames(vsd_reduced) <- c("HSPA9", "NKX3.1")
colnames(TCGA_reduced) <- c("HSPA9", "NKX3.1")

################################################################################
# SECTION 4: Correlation Analysis - Normalized Expression
################################################################################

# Perform Pearson correlation test on variance-stabilized data
res_vsd <- cor.test(vsd_reduced$HSPA9, vsd_reduced$NKX3.1, method = "pearson")

# Save correlation results
sink('../results/Co-expression/TCGA_correlation_normalized_statistics.txt')
cat("=== TCGA PRAD: NKX3.1 and HSPA9 Correlation (Normalized Expression) ===\n\n")
cat("Number of samples:", nrow(vsd_reduced), "\n\n")
cat("Pearson Correlation Test\n")
cat("------------------------\n\n")
print(res_vsd)
sink()

# Generate scatter plot with regression line (normalized data)
p_normalized <- ggscatter(
  vsd_reduced,
  x = "NKX3.1",
  y = "HSPA9",
  fullrange = TRUE,
  add = "reg.line",
  conf.int = FALSE,
  size = 0.5,
  cor.coef.coord = c(13.5, 14.75),
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
  ggplot2::annotate("text", x = 14.4, y = 14.6,
                    label = paste("n =", nrow(vsd_reduced)), size = 4) +
  xlab("NKX3.1 Expression (VST normalized)") +
  ylab("HSPA9 Expression (VST normalized)") +
  ggtitle("TCGA PRAD: NKX3.1 vs HSPA9 (Normalized Expression)")

# Save plot
ggsave("../results/Co-expression/TCGA_correlation_normalized.pdf",
       plot = p_normalized, height = 5, width = 5)

################################################################################
# SECTION 5: Correlation Analysis - Raw Counts
################################################################################

# Perform Pearson correlation test on raw counts
res_raw <- cor.test(TCGA_reduced$HSPA9, TCGA_reduced$NKX3.1, method = "pearson")

# Save correlation results
sink('../results/Co-expression/TCGA_correlation_raw_statistics.txt')
cat("=== TCGA PRAD: NKX3.1 and HSPA9 Correlation (Raw Counts) ===\n\n")
cat("Number of samples:", nrow(TCGA_reduced), "\n\n")
cat("Pearson Correlation Test\n")
cat("------------------------\n\n")
print(res_raw)
sink()

# Generate scatter plot with regression line (raw counts)
p_raw <- ggscatter(
  TCGA_reduced,
  x = "NKX3.1",
  y = "HSPA9",
  fullrange = TRUE,
  add = "reg.line",
  conf.int = FALSE,
  size = 0.5,
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
  ggplot2::annotate("text", x = 50000, y = 50000,
                    label = paste("n =", nrow(TCGA_reduced)), size = 4) +
  xlab("NKX3.1 Expression (Raw Counts)") +
  ylab("HSPA9 Expression (Raw Counts)") +
  ggtitle("TCGA PRAD: NKX3.1 vs HSPA9 (Raw Counts)")

# Save plot
ggsave("../results/Co-expression/TCGA_correlation_raw.pdf",
       plot = p_raw, height = 5, width = 5)

################################################################################
# SECTION 6: Save Summary Statistics
################################################################################

# Create summary data frame with correlation statistics
summary_df <- data.frame(
  Dataset = c("TCGA PRAD (Normalized)", "TCGA PRAD (Raw Counts)"),
  N_Samples = c(nrow(vsd_reduced), nrow(TCGA_reduced)),
  Pearson_R = c(res_vsd$estimate, res_raw$estimate),
  P_value = c(res_vsd$p.value, res_raw$p.value),
  CI_Lower = c(res_vsd$conf.int[1], res_raw$conf.int[1]),
  CI_Upper = c(res_vsd$conf.int[2], res_raw$conf.int[2])
)

write.csv(summary_df, "../results/Co-expression/TCGA_correlation_summary.csv",
          row.names = FALSE)

################################################################################
# Analysis Complete
################################################################################

cat("\nTCGA co-expression analysis completed successfully.\n")
cat("Results saved to: ../results/Co-expression/\n")
cat("\nCorrelation Summary:\n")
cat("  Normalized: R =", round(res_vsd$estimate, 3),
    ", p =", format(res_vsd$p.value, scientific = TRUE), "\n")
cat("  Raw Counts: R =", round(res_raw$estimate, 3),
    ", p =", format(res_raw$p.value, scientific = TRUE), "\n")
