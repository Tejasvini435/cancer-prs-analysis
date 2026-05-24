# Script: 01_gwas_qc_prs.R
# Project: Breast Cancer PRS Analysis
# Data: Michailidou et al. 2017 (Nature) GWAS summary statistics
# Author: Tejasvini Ramesh
# Date: May 2025
# Objective: To clean raw GWAS data, perform QC, Manhattan and QQ plots and export of clean SNP table 

library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)

filter <- dplyr::filter

# ---- SECTION 1: LOAD DATA ----------------------------------
gwas_raw <- fread("Projects/cancer-prs-analysis/data/raw/icogs_bcac_public_results_euro.txt.gz")

# Examining the data before performing any cleaning
dim(gwas_raw)        # how many rows and columns?
head(gwas_raw)       # what do the first 6 rows look like?
colnames(gwas_raw)   # what are the column names?

gwas_clean <- gwas_raw %>%
  rename(
    bp37 = position_b37,
    beta = bcac_icogs1_risk_beta,
    se   = bcac_icogs1_risk_se,
    p    = bcac_icogs1_risk_P1df,
    eaf  = bcac_icogs1_european_controls_eaf
  ) %>%

# QC filters
# Filter 1: Remove SNPs with missing values in key columns
filter(!is.na(beta), !is.na(se), !is.na(p), !is.na(eaf)) %>%
# Filter 2: Minor allele frequency > 1%
# Very rare variants are unreliable — small errors have big effects  
filter(eaf > 0.01 & eaf < 0.99) %>%
# Filter 3: Remove SNPs where SE is 0 or negative (data errors)
filter(se > 0) %>%
  # Filter 4: Remove strand-ambiguous SNPs (A/T and C/G — we can't tell which strand)
filter(!(a1 == "A" & a0 == "T") &
          !(a1 == "T" & a0 == "A") &
          !(a1 == "C" & a0 == "G") &
          !(a1 == "G" & a0 == "C")) %>%
  mutate(snp = paste0("chr", chr, ":", bp37))

print(paste0("SNPs after QC: ", nrow(gwas_clean)))
print(paste0("SNPs removed: ", nrow(gwas_raw) - nrow(gwas_clean)))

# ---- SECTION 3: GENOME-WIDE SIGNIFICANT HITS ---------------

# Genome-wide significance threshold is p < 5×10^-8
# This is a convention in GWAS — it corrects for testing millions of SNPs
gw_threshold <- 5e-8

gwsig <- gwas_clean %>% filter(p < gw_threshold)
print(paste0("Genome-wide significant SNPs: ", nrow(gwsig)))

# Look at the top 20 hits
top_hits <- gwas_clean %>%
  arrange(p) %>%
  head(20) %>%
  select(chr, bp37, snp, a1, beta, se, p, eaf)

print(top_hits)

# SECTION 4: MANHATTAN PLOT 
# Each dot is a SNP. X-axis = position in genome. Y-axis = -log10(p-value).
# The dots that poke above the red line are genome-wide significant.
# Manhattan needs a data frame with exactly these column names: SNP, CHR, BP, P (qqman package accepts only uppercase and letters, no numbers for the colnames)
manhattan_input <- gwas_clean %>%
  select(snp, chr, bp37, p) %>%
  rename(
    SNP = snp,
    CHR = chr,
    BP  = bp37,
    P   = p
  ) %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P)) %>%
  filter(CHR %in% 1:22) %>%
  mutate(
    CHR = as.integer(CHR),
    BP  = as.integer(BP)
  )

png("Projects/cancer-prs-analysis/results/manhattan_plot.png", width = 1400, height = 600, res = 120)
manhattan(
  manhattan_input,
  main      = "Breast Cancer GWAS - Michailidou et al. 2017 (Nature)",
  col       = c("#2E86AB", "#A23B72"),
  suggestiveline = -log10(1e-6),   # blue dotted line
  genomewideline = -log10(5e-8),   # red line
  cex       = 0.4,
  cex.axis  = 0.8
)
dev.off()
print(paste0("Manhattan plot saved."))

# SECTION 5: QQ PLOT
# A QQ plot checks whether your p-values are inflated.
# The dots should follow the diagonal line until the tail, where true hits deviate.
# If the whole line deviates upward, there's systematic bias which is bad.

png("Projects/cancer-prs-analysis/results/qq_plot.png", width = 600, height = 600, res = 120)
qq(manhattan_input$P,
   main = "QQ Plot — Breast Cancer GWAS (Michailidou 2017)")
dev.off()
print("QQ plot saved.")

# Section 6: Effect Size Distribution 
p1 <- ggplot(gwsig, aes(x = beta)) +
  geom_histogram(bins = 50, fill = "#2E86AB", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  labs(
    title    = "Effect size distribution of genome-wide significant SNPs",
    subtitle = paste0("n = ", nrow(gwsig), " SNPs, p < 5×10⁻⁸"),
    x        = "Effect size (log OR, beta)",
    y        = "Count",
    caption  = "Data: Michailidou et al. 2017, Nature"
  ) +
  theme_minimal(base_size = 13)

ggsave("Projects/cancer-prs-analysis/results/effect_size_distribution.png",
       p1, width = 8, height = 5, dpi = 150)

# SECTION 7: PRS THRESHOLD SENSITIVITY
# We test 5 different p-value cutoffs for including SNPs in the PRS.
# More SNPs = captures more signal but also more noise.

thresholds <- c(5e-8, 1e-6, 1e-4, 0.01, 0.05)
snp_counts <- sapply(thresholds, function(t) sum(gwas_clean$p < t, na.rm = TRUE))

threshold_df <- data.frame(
  threshold     = thresholds,
  log_threshold = -log10(thresholds),
  n_snps        = snp_counts
)

p2 <- ggplot(threshold_df, aes(x = factor(threshold), y = n_snps)) +
  geom_col(fill = "#A23B72", colour = "white") +
  geom_text(aes(label = format(n_snps, big.mark = ",")),
            vjust = -0.4, size = 4) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "PRS threshold sensitivity: SNP count by p-value cutoff",
    subtitle = "Breast Cancer GWAS - Michailidou et al. 2017",
    x        = "P-value threshold",
    y        = "Number of SNPs included in PRS",
    caption  = "More permissive thresholds capture polygenic signal but increase noise"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("Projects/cancer-prs-analysis/results/threshold_sensitivity.png", p2,
       width = 8, height = 5, dpi = 150)

# SECTION 8: CHROMOSOME SUMMARY

chrom_summary <- gwsig %>%
  group_by(chr) %>%
  summarise(
    n_snps    = n(),
    mean_beta = round(mean(beta), 4),
    mean_eaf  = round(mean(eaf), 3),
    min_p     = min(p)
  ) %>%
  arrange(as.integer(chr))

print(chrom_summary)
fwrite(chrom_summary, "Projects/cancer-prs-analysis/data/processed/chromosome_summary.csv")

# SECTION 9: EXPORT FOR PYTHON
# We export two clean files so the Python notebook can pick up where R left off.

# File 1: All QC-passed SNPs (for analysis in Python)
prs_export <- gwas_clean %>%
  filter(p < 0.05) %>%
  select(snp, chr, bp37, a1, a0, beta, se, p, eaf) %>%
  arrange(p)

fwrite(prs_export, "Projects/cancer-prs-analysis/data/processed/prs_snps_clean.csv")
print(paste0("Exported ", nrow(prs_export), " SNPs for Python analysis"))

# File 2: Threshold sensitivity table
fwrite(threshold_df, "Projects/cancer-prs-analysis/data/processed/threshold_sensitivity.csv")

print(paste0("R analysis complete. All outputs saved to results/ and data/processed/"))