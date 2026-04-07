01.run_gwas.sh        # Perform GWAS using GCTA with MLMA (Additive model GWAS main pipeline)
02.meta_GWAS.sh       # Conduct meta-analysis of GWAS results across multiple populations using METAL
03.QTL_define.R       # Define QTL regions from GWAS summary statistics, extract significant association peaks
04.finemapping.R      # Perform fine-mapping (using susieR) to identify potential causal variants within QTL regions
05.RNA_seq_processing.sh               # Automated pipeline for processing paired-end RNA-Seq data, performing: Quality control with fastp, Genome alignment using STAR, Gene quantification with featureCounts.
06.TPM_TMM_calculation_and_gene_filter.R                   # Performs TPM and TMM normalization, gene filtering.
07.peer_calculation.R              # Performs PEER factor analysis.
08.eqtl_mapping.sh                 # Performs cis/trans-eQTL analysis and eGene detection using TensorQTL.
09.coloc.R                         # Performs colocalization analysis between GWAS and eQTL datasets to identify shared genetic associations.

tensorqtl_parquet2txt.py           # Converts parquet files.
trans_eqtl.py                      # Performs trans-eQTL mapping using TensorQTL 
metal.txt             # METAL input configuration file specifying meta-analysis parameters and input GWAS summary files
