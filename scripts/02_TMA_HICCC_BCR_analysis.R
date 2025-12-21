################################################################################
# TMA HICCC - Biochemical Recurrence (BCR) Survival Analysis
################################################################################
#
# Project: NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer
#          Initiation
# PI: Cory Abate-Shen
# Publication: Cancer Discovery (2021)
# DOI: https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307
#
# Description:
# This script analyzes biochemical recurrence-free survival in the HICCC
# tissue microarray (TMA) cohort based on NKX3.1 subcellular localization
# (cytoplasmic, nuclear, and total expression).
#
################################################################################

# Load required libraries
library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

################################################################################
# SECTION 1: Helper Function - Read All Excel Sheets
################################################################################

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # Reads all sheets from an Excel file and returns a named list
  # Args:
  #   filename: path to Excel file
  #   tibble: if TRUE, returns tibbles; if FALSE, returns data.frames
  # Returns:
  #   Named list where names are sheet names and values are data frames

  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

################################################################################
# SECTION 2: Load TMA Data
################################################################################

mysheets <- read_excel_allsheets("../data/HICCC.xlsx")

# Create output directories
dir.create("../results/TMA-HICCC/BCR", recursive = TRUE, showWarnings = FALSE)

################################################################################
# SECTION 3: Cytoplasmic NKX3.1 BCR Analysis
################################################################################

# Prepare cytoplasmic NKX3.1 data
cyto_BCR <- mysheets$`Cytoplasmic NKX3.1 expression` %>%
  left_join(mysheets$`BCR-free Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(cyto_BCR)[2] <- "NKX3.1"

# Cox proportional hazards model
res.cox <- coxph(Surv(Months, BCR) ~ NKX3.1, data = cyto_BCR)

# Save statistical results
sink('../results/TMA-HICCC/BCR/cyto_BCR_statistics.txt')
cat("=== Cytoplasmic NKX3.1 - BCR Analysis ===\n\n")
cat("Cox Proportional Hazards Model\n")
cat("-------------------------------\n\n")
print(summary(res.cox))
sink()

# Fit survival curves
fit <- survfit(Surv(Months, BCR) ~ NKX3.1, data = cyto_BCR)

# Generate Kaplan-Meier plot with p-value
cyto_BCR_single <- ggsurvplot(
  fit,
  data = cyto_BCR,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

cyto_BCR_single$table <- cyto_BCR_single$table +
  labs(title = "Number of patients at risk")

# Add p-value annotation
cyto_BCR_single$plot <- cyto_BCR_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

# Save plot
pdf("../results/TMA-HICCC/BCR/cyto_BCR_KM_plot.pdf", width = 7, height = 9)
print(cyto_BCR_single, newpage = FALSE)
dev.off()

# Save summary statistics
cyto_summary <- cyto_BCR %>%
  group_by(NKX3.1) %>%
  summarize(
    n = n(),
    bcr_events = sum(BCR, na.rm = TRUE),
    .groups = 'drop'
  )
write.csv(cyto_summary, "../results/TMA-HICCC/BCR/cyto_BCR_summary.csv", row.names = FALSE)

################################################################################
# SECTION 4: Nuclear NKX3.1 BCR Analysis
################################################################################

# Prepare nuclear NKX3.1 data
nucl_BCR <- mysheets$`Nuclear NKX3.1 expression` %>%
  left_join(mysheets$`BCR-free Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(nucl_BCR)[2] <- "NKX3.1"

# Cox proportional hazards model
res.cox <- coxph(Surv(Months, BCR) ~ NKX3.1, data = nucl_BCR)

# Save statistical results
sink('../results/TMA-HICCC/BCR/nucl_BCR_statistics.txt')
cat("=== Nuclear NKX3.1 - BCR Analysis ===\n\n")
cat("Cox Proportional Hazards Model\n")
cat("-------------------------------\n\n")
print(summary(res.cox))
sink()

# Fit survival curves
fit <- survfit(Surv(Months, BCR) ~ NKX3.1, data = nucl_BCR)

# Generate Kaplan-Meier plot
nucl_BCR_single <- ggsurvplot(
  fit,
  data = nucl_BCR,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

nucl_BCR_single$table <- nucl_BCR_single$table +
  labs(title = "Number of patients at risk")

# Add p-value annotation
nucl_BCR_single$plot <- nucl_BCR_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

# Save plot
pdf("../results/TMA-HICCC/BCR/nucl_BCR_KM_plot.pdf", width = 7, height = 9)
print(nucl_BCR_single, newpage = FALSE)
dev.off()

# Save summary statistics
nucl_summary <- nucl_BCR %>%
  group_by(NKX3.1) %>%
  summarize(
    n = n(),
    bcr_events = sum(BCR, na.rm = TRUE),
    .groups = 'drop'
  )
write.csv(nucl_summary, "../results/TMA-HICCC/BCR/nucl_BCR_summary.csv", row.names = FALSE)

################################################################################
# SECTION 5: Total NKX3.1 BCR Analysis
################################################################################

# Prepare total NKX3.1 data
total_BCR <- mysheets$`Total NKX3.1 expression` %>%
  left_join(mysheets$`BCR-free Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(total_BCR)[2] <- "NKX3.1"

# Cox proportional hazards model
res.cox <- coxph(Surv(Months, BCR) ~ NKX3.1, data = total_BCR)

# Save statistical results
sink('../results/TMA-HICCC/BCR/total_BCR_statistics.txt')
cat("=== Total NKX3.1 - BCR Analysis ===\n\n")
cat("Cox Proportional Hazards Model\n")
cat("-------------------------------\n\n")
print(summary(res.cox))
sink()

# Fit survival curves
fit <- survfit(Surv(Months, BCR) ~ NKX3.1, data = total_BCR)

# Generate Kaplan-Meier plot
total_BCR_single <- ggsurvplot(
  fit,
  data = total_BCR,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("BCR-free estimated survival probability")) +
  xlab(c("Time (months)"))

total_BCR_single$table <- total_BCR_single$table +
  labs(title = "Number of patients at risk")

# Add p-value annotation
total_BCR_single$plot <- total_BCR_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

# Save plot
pdf("../results/TMA-HICCC/BCR/total_BCR_KM_plot.pdf", width = 7, height = 9)
print(total_BCR_single, newpage = FALSE)
dev.off()

# Save summary statistics
total_summary <- total_BCR %>%
  group_by(NKX3.1) %>%
  summarize(
    n = n(),
    bcr_events = sum(BCR, na.rm = TRUE),
    .groups = 'drop'
  )
write.csv(total_summary, "../results/TMA-HICCC/BCR/total_BCR_summary.csv", row.names = FALSE)

################################################################################
# Analysis Complete
################################################################################

cat("\nTMA-HICCC BCR survival analysis completed successfully.\n")
cat("Results saved to: ../results/TMA-HICCC/BCR/\n")
