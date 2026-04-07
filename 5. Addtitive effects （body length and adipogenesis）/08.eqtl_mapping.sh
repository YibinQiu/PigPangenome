#!/bin/bash
# ============================================================
# tensorqtl eQTL mapping 
# ============================================================

# --------------------------
# Section 1: Configuration Variables
# --------------------------

# Input files
GENO_PREFIX="samples_geno"               # PLINK genotype prefix
EXPR_BED="AF_TSS_TMM.bed.gz"            # Expression data (BED format)
COV_FILE="AF_cov.txt"                    # Covariates file

# Output prefixes
PERM_PREFIX="AF_perm"                    # Permutation results prefix
NOMINAL_PREFIX="AF_nominal"              # Nominal results prefix

# Chromosome range
CHROMOSOMES={1..18} 

# Script paths
PARQUET_CONVERTER="tensorqtl_parquet2txt.py"
EGENE_DETECTOR="tensorqtl_eGene_detection.R"
TRANS_ANALYZER="trans_eqtl.py"

# Output files
TRANS_RESULTS="trans_qtl_pv_filter.txt"  # trans-eQTL results
FDR_RESULTS="${PERM_PREFIX}.cis_qtl_adj_BH_fdr.txt.gz"  # FDR-adjusted results

# Log files
R_OUT_LOG="tensorqtl_eGene_detection.R.out"
R_ERR_LOG="tensorqtl_eGene_detection.R.err"

# --------------------------
# Section 2: Permutation cis-eQTL mapping
# --------------------------
echo "1. Running permutation cis-eQTL mapping..."
python3 -m tensorqtl \
    ${GENO_PREFIX} \
    ${EXPR_BED} \
    ${PERM_PREFIX} \
    --covariates ${COV_FILE} \
    --mode cis

# --------------------------
# Section 3: Nominal cis-eQTL mapping
# --------------------------
echo "2. Running nominal cis-eQTL mapping..."
python3 -m tensorqtl \
    ${GENO_PREFIX} \
    ${EXPR_BED} \
    ${NOMINAL_PREFIX} \
    --covariates ${COV_FILE} \
    --mode cis_nominal

# --------------------------
# Section 4: Parquet to Text Conversion
# --------------------------
echo "3. Converting parquet files to compressed text format..."
for CHR in "${CHROMOSOMES[@]}"; do
    INPUT_PARQUET="${NOMINAL_PREFIX}.cis_qtl_pairs.${CHR}.parquet"
    OUTPUT_TXT="${NOMINAL_PREFIX}.cis_qtl_pairs.${CHR}.txt.gz"
    
    echo "  Processing chromosome ${CHR}..."
    python3 ${PARQUET_CONVERTER} \
        "${INPUT_PARQUET}" \
        "${OUTPUT_TXT}"
done

# --------------------------
# Section 5: eGene Detection
# --------------------------
echo "4. Running eGene detection analysis..."
Rscript ${EGENE_DETECTOR} \
    "${PERM_PREFIX}.cis_qtl.txt.gz" \
    "${FDR_RESULTS}" \
    > "${R_OUT_LOG}" \
    2> "${R_ERR_LOG}"

# --------------------------
# Section 6: trans-eQTL Analysis
# --------------------------
echo "5. Running trans-eQTL analysis..."
python ${TRANS_ANALYZER} \
    ${GENO_PREFIX} \
    ${EXPR_BED} \
    ${COV_FILE} \
    ${TRANS_RESULTS}

echo "Process finished..."