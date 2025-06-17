#!/bin/bash
# usage:yhbatch 9.bench_ONT_vs_Chromosome_Scaffold_Contig.sbatch
# Job name:
#SBATCH --job-name=9.bench_ONT_vs_Chromosome_Scaffold_Contig.sbatch
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=rhenv
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=24
# Memory per node:
#SBATCH --output=9.bench_ONT_vs_Chromosome_Scaffold_Contig_%j.out
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C

module load anaconda3
module load bcftools
module load vcftools
source activate truvari
output=/BIGDATA2/scau_jyang_1/pangenome/projects/2023_PigGraphgenome_phase1/3.graph_genome_construct/scratch11
vcftools --gzvcf ${output}/2.ONT_merge.vcf.gz --indv 124 --recode --recode-INFO-all --stdout | awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} {$10="1/1"; print}' | bcftools sort | bgzip -f > ${output}/3.ONT_first_individual_modified.vcf.gz
tabix -p vcf ${output}/3.ONT_first_individual_modified.vcf.gz

vcftools --gzvcf ${output}/2.Chromosome_Scaffold_Contig_merge.vcf.gz --indv Bamei_ctg --recode --recode-INFO-all --stdout | awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} {$10="1/1"; print}' | bcftools sort | bgzip -f > ${output}/3.Chromosome_Scaffold_Contig_first_individual_modified.vcf.gz
tabix -p vcf ${output}/3.Chromosome_Scaffold_Contig_first_individual_modified.vcf.gz

truvari bench -b ${output}/3.ONT_first_individual_modified.vcf.gz -c ${output}/3.Chromosome_Scaffold_Contig_first_individual_modified.vcf.gz -O 0.0 -r 1000 -p 0.0 -P 0.3 -C 1000 -s 30 -S 10 -B 30 --sizemax 5000000 --pick multi --no-ref c -o ${output}/3.ONT_vs_Chromosome_Scaffold_Contig/

