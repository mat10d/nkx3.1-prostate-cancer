################################################################################
# TMA HICCC - Overall Survival Analysis
################################################################################
#
# Project: NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer
#          Initiation
# PI: Cory Abate-Shen
# Publication: Cancer Discovery (2021)
# DOI: https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307
#
# Description:
# This script analyzes overall survival in the HICCC tissue microarray (TMA)
# cohort based on NKX3.1 subcellular localization (cytoplasmic, nuclear, and
# total expression). It also performs multivariate analysis adjusting for
# clinical covariates (Gleason score, T-stage, age, SVI, PSA).
#
################################################################################

# Load required libraries
library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(cowplot)
library(grid)

################################################################################
# SECTION 1: Helper Function - Read All Excel Sheets
################################################################################

read_excel_allsheets <- function(filename, tibble = FALSE) {
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

# Create output directory
dir.create("../results/TMA-HICCC/Overall", recursive = TRUE, showWarnings = FALSE)

################################################################################
# SECTION 3: Cytoplasmic NKX3.1 Overall Survival Analysis
################################################################################

# Prepare cytoplasmic NKX3.1 data
cyto_overall <- mysheets$`Cytoplasmic NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(cyto_overall)[2] <- "NKX3.1"

# Cox proportional hazards model
res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1, data = cyto_overall)

# Save statistical results
sink('../results/TMA-HICCC/Overall/cyto_overall_statistics.txt')
cat("=== Cytoplasmic NKX3.1 - Overall Survival Analysis ===\n\n")
print(summary(res.cox))
sink()

# Fit survival curves
fit <- survfit(Surv(Months, Dead) ~ NKX3.1, data = cyto_overall)

# Generate Kaplan-Meier plot
cyto_overall_single <- ggsurvplot(
  fit,
  data = cyto_overall,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("Overall Survival Probability")) +
  xlab(c("Time (months)"))

cyto_overall_single$table <- cyto_overall_single$table +
  labs(title = "Number of patients at risk")

cyto_overall_single$plot <- cyto_overall_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

pdf("../results/TMA-HICCC/Overall/cyto_overall_KM_plot.pdf", width = 7, height = 9)
print(cyto_overall_single, newpage = FALSE)
dev.off()

# Save summary statistics
cyto_summary <- cyto_overall %>%
  group_by(NKX3.1) %>%
  summarize(n = n(), deaths = sum(Dead, na.rm = TRUE), .groups = 'drop')
write.csv(cyto_summary, "../results/TMA-HICCC/Overall/cyto_overall_summary.csv", row.names = FALSE)

################################################################################
# SECTION 4: Nuclear NKX3.1 Overall Survival Analysis
################################################################################

nucl_overall <- mysheets$`Nuclear NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(nucl_overall)[2] <- "NKX3.1"

res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1, data = nucl_overall)

sink('../results/TMA-HICCC/Overall/nucl_overall_statistics.txt')
cat("=== Nuclear NKX3.1 - Overall Survival Analysis ===\n\n")
print(summary(res.cox))
sink()

fit <- survfit(Surv(Months, Dead) ~ NKX3.1, data = nucl_overall)

nucl_overall_single <- ggsurvplot(
  fit,
  data = nucl_overall,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("Overall Survival Probability")) +
  xlab(c("Time (months)"))

nucl_overall_single$table <- nucl_overall_single$table +
  labs(title = "Number of patients at risk")

nucl_overall_single$plot <- nucl_overall_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

pdf("../results/TMA-HICCC/Overall/nucl_overall_KM_plot.pdf", width = 7, height = 9)
print(nucl_overall_single, newpage = FALSE)
dev.off()

nucl_summary <- nucl_overall %>%
  group_by(NKX3.1) %>%
  summarize(n = n(), deaths = sum(Dead, na.rm = TRUE), .groups = 'drop')
write.csv(nucl_summary, "../results/TMA-HICCC/Overall/nucl_overall_summary.csv", row.names = FALSE)

################################################################################
# SECTION 5: Total NKX3.1 Overall Survival Analysis
################################################################################

total_overall <- mysheets$`Total NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1)

names(total_overall)[2] <- "NKX3.1"

res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1, data = total_overall)

sink('../results/TMA-HICCC/Overall/total_overall_statistics.txt')
cat("=== Total NKX3.1 - Overall Survival Analysis ===\n\n")
print(summary(res.cox))
sink()

fit <- survfit(Surv(Months, Dead) ~ NKX3.1, data = total_overall)

total_overall_single <- ggsurvplot(
  fit,
  data = total_overall,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "red"),
  legend.labs = c("Low NKX3.1", "High NKX3.1"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("Overall Survival Probability")) +
  xlab(c("Time (months)"))

total_overall_single$table <- total_overall_single$table +
  labs(title = "Number of patients at risk")

total_overall_single$plot <- total_overall_single$plot +
  ggplot2::annotate("text", x = 50, y = 0.30,
                    label = paste("Log Rank p-value =",
                                  round(as.numeric(surv_pvalue(fit)[2]), digits = 4)))

pdf("../results/TMA-HICCC/Overall/total_overall_KM_plot.pdf", width = 7, height = 9)
print(total_overall_single, newpage = FALSE)
dev.off()

total_summary <- total_overall %>%
  group_by(NKX3.1) %>%
  summarize(n = n(), deaths = sum(Dead, na.rm = TRUE), .groups = 'drop')
write.csv(total_summary, "../results/TMA-HICCC/Overall/total_overall_summary.csv", row.names = FALSE)

################################################################################
# SECTION 6: Combined Nuclear and Cytoplasmic Analysis
################################################################################

nuc_cyto_overall <- mysheets$`Nuclear NKX3.1 expression` %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1) %>%
  dplyr::rename(Nucl_NKX3.1 = `NKX3.1 High`) %>%
  dplyr::select(-`NKX3.1 Low`) %>%
  left_join(mysheets$`Cytoplasmic NKX3.1 expression`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1) %>%
  dplyr::rename(Cyto_NKX3.1 = `NKX3.1 High`) %>%
  dplyr::select(-`NKX3.1 Low`) %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #")

fit <- survfit(Surv(Months, Dead) ~ Nucl_NKX3.1 + Cyto_NKX3.1, data = nuc_cyto_overall)

nuc_cyto_single <- ggsurvplot(
  fit,
  data = nuc_cyto_overall,
  conf.int = FALSE,
  risk.table = TRUE,
  palette = c("blue", "darkturquoise", "salmon", "red"),
  legend.labs = c("low nucl, low cyto",
                  "low nucl, high cyto",
                  "high nucl, low cyto",
                  "high nucl, high cyto"),
  pval = FALSE,
  ggtheme = theme_gray()
) +
  ylab(c("Overall Survival Probability")) +
  xlab(c("Time (months)"))

nuc_cyto_single$table <- nuc_cyto_single$table +
  labs(title = "Number of patients at risk") +
  ylab("Strata")

pdf("../results/TMA-HICCC/Overall/nuc_cyto_overall_combined_KM_plot.pdf", width = 8, height = 9)
print(nuc_cyto_single, newpage = FALSE)
dev.off()

res <- pairwise_survdiff(Surv(Months, Dead) ~ Nucl_NKX3.1 + Cyto_NKX3.1,
                         data = nuc_cyto_overall, p.adjust.method = "none")
res.cox <- coxph(Surv(Months, Dead) ~ Nucl_NKX3.1 + Cyto_NKX3.1, data = nuc_cyto_overall)

sink('../results/TMA-HICCC/Overall/nuc_cyto_overall_statistics.txt')
cat("=== Combined Nuclear and Cytoplasmic NKX3.1 - Overall Survival ===\n\n")
print(res)
cat("\n\n")
print(summary(res.cox))
sink()

################################################################################
# SECTION 7: Clinical Data Processing
################################################################################

# Process Gleason scores
mysheets$`Gleason Score`[is.na(mysheets$`Gleason Score`)] <- "N"

for(i in 1:nrow(mysheets$`Gleason Score`)){
  if(mysheets$`Gleason Score`$five[i] == "Y") mysheets$`Gleason Score`$g_score[i] <- 5
  if(mysheets$`Gleason Score`$six[i] == "Y") mysheets$`Gleason Score`$g_score[i] <- 6
  if(mysheets$`Gleason Score`$seven[i] == "Y") mysheets$`Gleason Score`$g_score[i] <- 7
  if(mysheets$`Gleason Score`$eight[i] == "Y") mysheets$`Gleason Score`$g_score[i] <- 8
  if(mysheets$`Gleason Score`$nine[i] == "Y") mysheets$`Gleason Score`$g_score[i] <- 9
}

# Process T-stage
mysheets$`T-stage`[is.na(mysheets$`T-stage`)] <- "N"

for(i in 1:nrow(mysheets$`T-stage`)){
  if(mysheets$`T-stage`$T2a[i] == "Y") mysheets$`T-stage`$t_score[i] <- 2
  if(mysheets$`T-stage`$T2c[i] == "Y") mysheets$`T-stage`$t_score[i] <- 2
  if(mysheets$`T-stage`$T3a[i] == "Y") mysheets$`T-stage`$t_score[i] <- 3
  if(mysheets$`T-stage`$T3b[i] == "Y") mysheets$`T-stage`$t_score[i] <- 3
  if(mysheets$`T-stage`$T3c[i] == "Y") mysheets$`T-stage`$t_score[i] <- 3
  if(mysheets$`T-stage`$T4[i] == "Y") mysheets$`T-stage`$t_score[i] <- 4
}

################################################################################
# SECTION 8: Multivariate Analysis - Cytoplasmic NKX3.1
################################################################################

# Univariate analysis helper function
perform_univariate_analysis <- function(data, covariates) {
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(Months, Dead) ~', x)))

  univ_models <- lapply(univ_formulas, function(x) coxph(x, data = data))

  univ_results <- lapply(univ_models, function(x) {
    x <- summary(x)
    p.value <- signif(x$wald["pvalue"], digits = 2)
    wald.test <- signif(x$wald["test"], digits = 2)
    beta <- signif(x$coef[1], digits = 2)
    HR <- signif(x$coef[2], digits = 2)
    HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
    HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
    res <- c(beta, HR, wald.test, p.value)
    names(res) <- c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
    return(res)
  })

  return(t(as.data.frame(univ_results, check.names = FALSE)))
}

# Cytoplasmic - Multivariate with full covariates
cyto_overall_multi <- mysheets$`Cytoplasmic NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  left_join(mysheets$`Gleason Score`, by = "TMA #") %>%
  left_join(mysheets$`T-stage`, by = "TMA #") %>%
  left_join(mysheets$`Age at RP`, by = "TMA #") %>%
  left_join(mysheets$SVI, by = "TMA #") %>%
  left_join(mysheets$`PSA levels`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1) %>%
  dplyr::rename(NKX3.1 = `NKX3.1 High`) %>%
  dplyr::select(NKX3.1, Months, Dead, g_score, t_score, Years, SVI, PSA) %>%
  dplyr::filter(complete.cases(.))

# Multivariate Cox model
res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1 + g_score + t_score + Years + SVI + PSA,
                 data = cyto_overall_multi)

sink('../results/TMA-HICCC/Overall/cyto_overall_multivariate_statistics.txt')
cat("=== Cytoplasmic NKX3.1 - Multivariate Analysis ===\n\n")
print(summary(res.cox))
sink()

# Univariate analysis for each covariate
covariates <- c("NKX3.1", "g_score", "t_score", "Years", "SVI", "PSA")
res <- perform_univariate_analysis(cyto_overall_multi, covariates)
write.csv(as.data.frame(res), "../results/TMA-HICCC/Overall/cyto_overall_univariate_results.csv")

################################################################################
# SECTION 9: Multivariate Analysis - Nuclear NKX3.1
################################################################################

nucl_overall_multi <- mysheets$`Nuclear NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  left_join(mysheets$`Gleason Score`, by = "TMA #") %>%
  left_join(mysheets$`T-stage`, by = "TMA #") %>%
  left_join(mysheets$`Age at RP`, by = "TMA #") %>%
  left_join(mysheets$SVI, by = "TMA #") %>%
  left_join(mysheets$`PSA levels`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1) %>%
  dplyr::rename(NKX3.1 = `NKX3.1 High`) %>%
  dplyr::select(NKX3.1, Months, Dead, g_score, t_score, Years, SVI, PSA) %>%
  dplyr::filter(complete.cases(.))

res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1 + g_score + t_score + Years + SVI + PSA,
                 data = nucl_overall_multi)

sink('../results/TMA-HICCC/Overall/nucl_overall_multivariate_statistics.txt')
cat("=== Nuclear NKX3.1 - Multivariate Analysis ===\n\n")
print(summary(res.cox))
sink()

res <- perform_univariate_analysis(nucl_overall_multi, covariates)
write.csv(as.data.frame(res), "../results/TMA-HICCC/Overall/nucl_overall_univariate_results.csv")

################################################################################
# SECTION 10: Multivariate Analysis - Total NKX3.1
################################################################################

total_overall_multi <- mysheets$`Total NKX3.1 expression` %>%
  left_join(mysheets$`Overall Survival`, by = "TMA #") %>%
  left_join(mysheets$`Gleason Score`, by = "TMA #") %>%
  left_join(mysheets$`T-stage`, by = "TMA #") %>%
  left_join(mysheets$`Age at RP`, by = "TMA #") %>%
  left_join(mysheets$SVI, by = "TMA #") %>%
  left_join(mysheets$`PSA levels`, by = "TMA #") %>%
  mutate_at(vars(`NKX3.1 High`), recode, "Y" = 2) %>%
  mutate_at(vars(`NKX3.1 High`), replace_na, 1) %>%
  dplyr::rename(NKX3.1 = `NKX3.1 High`) %>%
  dplyr::select(NKX3.1, Months, Dead, g_score, t_score, Years, SVI, PSA) %>%
  dplyr::filter(complete.cases(.))

res.cox <- coxph(Surv(Months, Dead) ~ NKX3.1 + g_score + t_score + Years + SVI + PSA,
                 data = total_overall_multi)

sink('../results/TMA-HICCC/Overall/total_overall_multivariate_statistics.txt')
cat("=== Total NKX3.1 - Multivariate Analysis ===\n\n")
print(summary(res.cox))
sink()

res <- perform_univariate_analysis(total_overall_multi, covariates)
write.csv(as.data.frame(res), "../results/TMA-HICCC/Overall/total_overall_univariate_results.csv")

################################################################################
# Analysis Complete
################################################################################

cat("\nTMA-HICCC Overall Survival analysis completed successfully.\n")
cat("Results saved to: ../results/TMA-HICCC/Overall/\n")
