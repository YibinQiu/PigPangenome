#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Argument parsing
snp_ids_file <- args[1]
ld_file <- args[2]
sumstats_file <- args[3]
sample_size <- 1496
output_prefix <- args[4]

# Load required package
library(susieR)

# Step 1: Read SNP IDs
snp_ids <- scan(snp_ids_file, what = character())

# Step 2: Read LD matrix
R <- as.matrix(read.table(ld_file))
colnames(R) <- rownames(R) <- snp_ids

# Step 3: Read GWAS summary statistics
sumstats <- read.table(sumstats_file, header = TRUE)

# Step 4: Match order 
if (!all(snp_ids %in% sumstats[,2])) {
  stop("Some SNPs in the LD matrix are not found in the GWAS summary. Check SNP ID formats.")
}
sumstats <- sumstats[match(snp_ids, sumstats[,2]), ]
z <- as.numeric(sumstats[,3])

snp_ordered <- sumstats[,2]

# Step 5: Run SuSiE using summary statistics and LD matrix
res <- susie_rss(z = z, R = R, L = 10, coverage = 0.95, n = sample_size)

# Step 6: Output credible sets (each set saved in a single line)
cs <- susie_get_cs(res)
cs_snps <- lapply(cs$cs, function(set) snp_ordered[set])
lines <- sapply(cs_snps, function(x) paste(x, collapse = "\t"))
writeLines(lines, paste0(output_prefix, "_credible_sets.txt"))

