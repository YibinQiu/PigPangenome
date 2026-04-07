# ==============================================================================
# Colocalization Analysis
# ==============================================================================

# --------------------------
# Section 1: Load Required Packages
# --------------------------
library(dplyr)
library(coloc)

# --------------------------
# Section 2: Handle Command Line Arguments
# --------------------------
args <- commandArgs(trailingOnly = TRUE)

# Define input/output parameters
GWAS_PVAL_FILE    <- args[1]  # GWAS p-values file
GWAS_SAMPLE_SIZE  <- args[2]  # GWAS sample size
EQTL_PVAL_FILE    <- args[3]  # eQTL p-values file 
EQTL_SAMPLE_SIZE  <- args[4]  # eQTL sample size
OUTPUT_FILE       <- args[5]  # Output file

# --------------------------
# Section 3: Data Preparation
# --------------------------

# Read and process GWAS data
gwas_data <- read.table(GWAS_PVAL_FILE, header = TRUE)
colnames(gwas_data) <- c("chr", "variant_id", "pos", "allele1", "allele0", 
                        "maf", "beta", "se", "pval_nominal", "pve", "variant_type")

gwas$varbeta <- gwas$se^2
gwas$N <- input_2
gwas$N <- as.numeric(gwas$N)
  
# Read and process eQTL data
eqtl_data <- read.table(EQTL_PVAL_FILE, header = TRUE)
colnames(eqtl_data) <- c("gene_id", "variant_id", "tss_distance", "maf", "ma_samples",
                        "ma_count", "pval_nominal", "slope", "slope_se", "variant_type")

eqtl$varbeta <- eqtl$slope_se^2
eqtl$N <- input_4
eqtl$N <- as.numeric(eqtl$N)

# Merge datasets
merged_data <- inner_join(eqtl_data, gwas_data, 
                         by = "variant_id", 
                         suffix = c("_eqtl", "_gwas"))

cat("Number of variants after merging:", nrow(merged_data), "\n")

# --------------------------
# Section 4: Colocalization Analysis
# --------------------------

# Run coloc analysis
coloc_results <- coloc.abf(
  dataset1 = list(
    snp = merged_data$variant_id,
    pvalues = merged_data$pval_nominal_gwas,
    N = merged_data$N_gwas,
    type = "quant",
    beta = merged_data$beta,
    varbeta = merged_data$varbeta_gwas,
    MAF = merged_data$maf_gwas
  ),
  dataset2 = list(
    snp = merged_data$variant_id,
    pvalues = merged_data$pval_nominal_eqtl,
    N = merged_data$N_eqtl,
    type = "quant",
    beta = merged_data$slope,
    varbeta = merged_data$varbeta_eqtl,
    MAF = merged_data$maf_eqtl
  ),
  p1 = 1/(nrow(merged_data) + 1),
  p2 = 1/(nrow(merged_data) + 1),
  p12 = 1e-5
)

# Filter results for significant colocalizations (PP.H4 > 0.05)
significant_results <- coloc_results$results %>%
  filter(SNP.PP.H4 > 0.05)

# --------------------------
# Section 5: Output Results
# --------------------------

# Write significant results
write.table(
  significant_results,
  file = OUTPUT_FILE,
  sep = " ",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("Colocalization analysis completed. Results saved to", OUTPUT_FILE, "\n")