00.replace_genotype.py  # Generate genotype files under dominant/recessive models by replacing heterozygotes (0/1)
01.run_gwas.sh          # Run GWAS analysis (dominant/recessive models) using GCTA with MLMA
02.meta_GWAS.sh         # Conduct meta-analysis of GWAS results across multiple populations using METAL
03.QTL_define.R         # Define significant QTL regions from meta-GWAS results
metal.txt               # Configuration file for METAL specifying input files and meta-analysis parameters
04.RNA_seq_processing.sh               # Automated pipeline for processing paired-end RNA-Seq data, performing: Quality control with fastp, Genome alignment using STAR, Gene quantification with featureCounts.
05.TPM_TMM_calculation_and_gene_filter.R                   # Performs TPM and TMM normalization, gene filtering.
06.peer_calculation.R              # Performs PEER factor analysis.
07.eqtl_mapping.sh                 # Performs cis/trans-eQTL analysis and eGene detection using TensorQTL.
08.coloc.R                         # Performs colocalization analysis between GWAS and eQTL datasets to identify shared genetic associations.
tensorqtl_parquet2txt.py           # Converts parquet files.
trans_eqtl.py                      # Performs trans-eQTL mapping using TensorQTL 