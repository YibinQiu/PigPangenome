# ============================================================
# TensorQTL trans-eQTL Mapping 
# ============================================================

import os
import re
import pandas as pd
import torch
import tensorqtl
import sys
from tensorqtl import pgen, cis, trans, post, genotypeio

# --------------------------
# Section 1: Setup and Configuration
# --------------------------

# GPU if available, otherwise CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")

# --------------------------
# Section 2: Parameter Handling
# --------------------------

# Get command line arguments
plink_prefix_path = sys.argv[1]  # Path to PLINK genotype files prefix
expression_bed = sys.argv[2]     # Path to expression BED file
covariates_file = sys.argv[3]    # Path to covariates file
output_prefix = sys.argv[4]      # Output file prefix

# Print parameters for verification
print(f"Genotype path: {plink_prefix_path}")
print(f"Expression BED: {expression_bed}")
print(f"Covariates file: {covariates_file}")
print(f"Output prefix: {output_prefix}")

# --------------------------
# Section 3: Data Loading
# --------------------------

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

covariates_df = pd.read_csv(covariates_file, 
                           sep='\t', 
                           engine='python', 
                           index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)

genotype_df = pr.load_genotypes()

variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# --------------------------
# Section 4: Trans-eQTL Mapping
# --------------------------

#run trans-eQTL mapping
#to limit output size, only associations with p-value <= 5e-04 are returned
trans_df = trans.map_trans(genotype_df, 
                          phenotype_df, 
                          covariates_df,
                          return_sparse=True, 
                          pval_threshold=5e-04, 
                          batch_size=20000, 
                          maf_threshold=0.05)

# --------------------------
# Section 5: Post-processing
# --------------------------

# Remove cis-associations (within 1Mb window of gene)
trans_df = trans.filter_cis(trans_df, 
                           phenotype_pos_df, 
                           variant_df, 
                           window=1000000)

# --------------------------
# Section 6: Output Results
# --------------------------

# Save trans-eQTL results
trans_df.to_csv(output_prefix, 
               sep='\t', 
               encoding='utf-8', 
               index=False, 
               header=True)