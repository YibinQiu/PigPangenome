#!/bin/bash

# ---------------- Parameters ----------------
if [ $# -lt 1 ]; then
  echo "Usage: sbatch 1.run_gwas.sh <population>"
  exit 1
fi

pop=$1
threads=40

# Paths
project=$HOME/pangenome/project/GWAS/dominance
genofile_path=$project/data/${pop}_merge

raw_dir=$project/${pop}_merge/raw_data
grm_dir=$project/${pop}_merge/GRM
pca_dir=$project/${pop}_merge/PCA
gwas_dir=$project/${pop}_merge/GWAS
mkdir -p $raw_dir $grm_dir $pca_dir $gwas_dir logs

# ---------------- Step 1: Convert genotype data to PLINK binary format ----------------
echo "[Step 1] Converting genotype data to PLINK binary format"
plink \
  --bfile ${genofile_path} \
  --keep ${raw_dir}/${pop}_samples.txt \
  --make-bed \
  --out ${raw_dir}/${pop}_merge \
  --allow-extra-chr \
  --chr 1-18 \
  --keep-allele-order \
  --threads ${threads}

# ---------------- Step 2: Generate the Genetic Relationship Matrix (GRM) ----------------
echo "[Step 2] Generating GRM"
gcta \
  --bfile ${raw_dir}/${pop}_merge \
  --make-grm \
  --autosome \
  --thread-num ${threads} \
  --out ${grm_dir}/${pop}_merge

# ---------------- Step 3: Calculate Principal Components (PCA) ----------------
echo "[Step 3] Calculating PCA"
gcta \
  --grm ${grm_dir}/${pop}_merge \
  --pca 5 \
  --thread-num ${threads} \
  --out ${pca_dir}/${pop}_merge

# ---------------- Step 4: Perform GWAS using the MLMA model ----------------
echo "[Step 4] Running GWAS with MLMA model"
gcta \
  --mlma \
  --bfile ${raw_dir}/${pop}_merge \
  --grm ${grm_dir}/${pop}_merge \
  --pheno ${raw_dir}/${pop}_pheno.txt \
  --mpheno 1 \
  --covar ${raw_dir}/${pop}_covar.txt \
  --qcovar ${raw_dir}/${pop}_qcovar.txt \
  --out ${gwas_dir}/${pop}_GWAS \
  --thread-num ${threads}

echo "[Done] GWAS completed for population: $pop"
