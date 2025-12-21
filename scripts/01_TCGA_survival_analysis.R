################################################################################
# TCGA Survival Analysis - NKX3.1 and HSPA9
################################################################################
#
# Project: NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer
#          Initiation
# PI: Cory Abate-Shen
# Publication: Cancer Discovery (2021)
# DOI: https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307
#
# Description:
# This script analyzes the relationship between HSPA9 and NKX3.1 expression
# and biochemical recurrence (BCR) in TCGA prostate cancer samples.
# It performs survival analysis using Cox proportional hazards models and
# generates Kaplan-Meier plots.
#
################################################################################

# Load required libraries
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(DESeq2)
library(cowplot)
library(grid)

################################################################################
# SECTION 1: Data Loading and Preparation
################################################################################

# Load TCGA processed data
TCGA <- readRDS("../data/TCGA-processed-data.rds")

# Extract survival-related variables
days_to_last_follow_up <- TCGA$recount_colData$gdc_cases.diagnoses.days_to_last_follow_up
time_to_recurrence <- TCGA$recount_colData$xml_days_to_first_biochemical_recurrence
recurrence <- TCGA$recount_colData$xml_biochemical_recurrence

# Create survival data frame with relevant columns
TCGA_survival <- TCGA$recount_colData %>%
  dplyr::select(
    gdc_cases.samples.submitter_id,
    gdc_cases.diagnoses.days_to_last_follow_up,
    xml_days_to_first_biochemical_recurrence,
    xml_biochemical_recurrence
  )

# Impute missing recurrence time with last follow-up time
for(i in 1:nrow(TCGA_survival)){
  if(is.na(TCGA_survival$xml_days_to_first_biochemical_recurrence[i])){
    TCGA_survival$xml_days_to_first_biochemical_recurrence[i] <-
      TCGA_survival$gdc_cases.diagnoses.days_to_last_follow_up[i]
  }
}

# Remove samples with missing data
TCGA_survival <- TCGA_survival[complete.cases(TCGA_survival), ]

################################################################################
# SECTION 2: Gene Expression Processing
################################################################################

# Filter counts to include only samples with survival data
counts_reduced <- TCGA$counts[, colnames(TCGA$counts) %in%
                                 TCGA_survival$gdc_cases.samples.submitter_id]

# Variance stabilizing transformation
vsd <- vst(TCGA$counts, blind = FALSE)

# Extract HSPA9 (ENSG00000167034) and NKX3.1 (ENSG00000113013) expression
vsd_reduced <- as.data.frame(t(vsd[rownames(vsd) %in%
                                     c("ENSG00000167034", "ENSG00000113013"), ]))
colnames(vsd_reduced) <- c("HSPA9", "NKX3.1")

################################################################################
# SECTION 3: Gene Expression Stratification
################################################################################

# Categorize samples based on gene expression levels
# HSPA9: median split (high vs low)
# HSPA9_extreme: quartile split (top 25% vs bottom 25%)
# NKX3.1: quartile split (top 25% vs bottom 25%)

for(i in 1:nrow(TCGA_survival)){
  sample <- TCGA_survival$gdc_cases.samples.submitter_id[i]

  # HSPA9 median split (2 = high, 1 = low)
  TCGA_survival$HSPA9[i] <- ifelse(
    vsd_reduced[sample, "HSPA9"] > quantile(vsd_reduced[, "HSPA9"])["50%"],
    2, 1
  )

  # HSPA9 extreme quartiles (2 = top 25%, 1 = bottom 25%, NA = middle 50%)
  TCGA_survival$HSPA9_extreme[i] <- ifelse(
    vsd_reduced[sample, "HSPA9"] > quantile(vsd_reduced[, "HSPA9"])["75%"],
    2,
    ifelse(vsd_reduced[sample, "HSPA9"] < quantile(vsd_reduced[, "HSPA9"])["25%"],
           1, NA)
  )

  # NKX3.1 quartile split (2 = top 25%, 1 = bottom 25%, NA = middle 50%)
  TCGA_survival$NKX3.1[i] <- ifelse(
    vsd_reduced[sample, "NKX3.1"] > quantile(vsd_reduced[, "NKX3.1"])["75%"],
    2,
    ifelse(vsd_reduced[sample, "NKX3.1"] < quantile(vsd_reduced[, "NKX3.1"])["25%"],
           1, NA)
  )
}

################################################################################
# SECTION 4: Combined HSPA9 and NKX3.1 Survival Analysis
################################################################################

# Prepare data for combined analysis
bcr_TCGA_survival <- TCGA_survival[complete.cases(TCGA_survival), ]
bcr_TCGA_survival$BCR <- ifelse(bcr_TCGA_survival$xml_biochemical_recurrence == "YES", 2, 1)

# Create output directory
dir.create("../results/TCGA", recursive = TRUE, showWarnings = FALSE)

# Fit survival curves for combined HSPA9 and NKX3.1 groups
fit <- survfit(Surv(xml_days_to_first_biochemical_recurrence, BCR) ~
                 HSPA9 + NKX3.1, data = bcr_TCGA_survival)

# Generate Kaplan-Meier plot (without p-value)
bcr_TCGA_single <- ggsurvplot(
  fit,
  data = bcr_TCGA_survival,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "darkturquoise", "salmon", "red"),
  legend.labs = c("low HSPA9, low NKX3.1",
                  "low HSPA9, high NKX3.1",
                  "high HSPA9, low NKX3.1",
                  "high HSPA9, high NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("BCR Survival Probability")) +
  xlab(c("Time (days)"))

bcr_TCGA_single$table <- bcr_TCGA_single$table +
  labs(title = "Number of patients at risk") +
  ylab("Strata")

# Save plot
pdf("../results/TCGA/HSPA9_NKX3.1_bcr_combined.pdf", width = 9, height = 9)
print(bcr_TCGA_single, newpage = FALSE)
dev.off()

# Perform pairwise survival comparison and Cox regression
res <- pairwise_survdiff(
  Surv(xml_days_to_first_biochemical_recurrence, BCR) ~ HSPA9 + NKX3.1,
  data = bcr_TCGA_survival,
  p.adjust.method = "none"
)

res.cox <- coxph(
  Surv(xml_days_to_first_biochemical_recurrence, BCR) ~ HSPA9 + NKX3.1,
  data = bcr_TCGA_survival
)

# Save statistical results
sink('../results/TCGA/HSPA9_NKX3.1_bcr_combined_statistics.txt')
cat("=== Pairwise Survival Comparison ===\n\n")
print(res)
cat("\n\n=== Cox Proportional Hazards Model ===\n\n")
print(summary(res.cox))
sink()

# Save summary statistics as CSV
bcr_summary <- bcr_TCGA_survival %>%
  group_by(HSPA9, NKX3.1) %>%
  summarize(
    n = n(),
    bcr_events = sum(BCR == 2),
    .groups = 'drop'
  )
write.csv(bcr_summary, "../results/TCGA/HSPA9_NKX3.1_sample_summary.csv", row.names = FALSE)

################################################################################
# SECTION 5: HSPA9-Only Survival Analysis (Extreme Quartiles)
################################################################################

# Prepare data for HSPA9 extreme quartiles analysis
bcr_TCGA_survival <- TCGA_survival
bcr_TCGA_survival$BCR <- ifelse(bcr_TCGA_survival$xml_biochemical_recurrence == "YES", 2, 1)

# Fit survival curves for HSPA9 extreme groups
fit <- survfit(Surv(xml_days_to_first_biochemical_recurrence, BCR) ~
                 HSPA9_extreme, data = bcr_TCGA_survival)

# Generate Kaplan-Meier plot
bcr_TCGA_single <- ggsurvplot(
  fit,
  data = bcr_TCGA_survival,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("low HSPA9", "high HSPA9"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("BCR Survival Probability")) +
  xlab(c("Time (days)"))

bcr_TCGA_single$table <- bcr_TCGA_single$table +
  labs(title = "Number of patients at risk") +
  ylab("Strata")

# Add p-value annotation
bcr_TCGA_single$plot <- bcr_TCGA_single$plot +
  ggplot2::annotate("text", x = 2000, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

# Save plot
pdf("../results/TCGA/HSPA9_extreme_bcr.pdf", width = 8, height = 9)
print(bcr_TCGA_single, newpage = FALSE)
dev.off()

# Perform Cox regression
res.cox <- coxph(
  Surv(xml_days_to_first_biochemical_recurrence, BCR) ~ HSPA9_extreme,
  data = bcr_TCGA_survival
)

# Save statistical results
sink('../results/TCGA/HSPA9_extreme_bcr_statistics.txt')
cat("=== Cox Proportional Hazards Model (HSPA9 Extreme Quartiles) ===\n\n")
print(summary(res.cox))
sink()

# Save summary statistics as CSV
hspa9_summary <- bcr_TCGA_survival %>%
  filter(!is.na(HSPA9_extreme)) %>%
  group_by(HSPA9_extreme) %>%
  summarize(
    n = n(),
    bcr_events = sum(BCR == 2, na.rm = TRUE),
    .groups = 'drop'
  )
write.csv(hspa9_summary, "../results/TCGA/HSPA9_extreme_sample_summary.csv", row.names = FALSE)

################################################################################
# Analysis Complete
################################################################################

cat("\nTCGA survival analysis completed successfully.\n")
cat("Results saved to: ../results/TCGA/\n")
