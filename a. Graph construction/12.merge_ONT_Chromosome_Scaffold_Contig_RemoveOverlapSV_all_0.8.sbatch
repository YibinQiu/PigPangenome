#!/bin/bash
# usage:yhbatch 12.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV_all_0.8.sbatch
# Job name:
#SBATCH --job-name=12.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV_all_0.8.sbatch
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
#SBATCH --output=12.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV_all_0.8_%j.out
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C

module load anaconda3
module load bcftools
module load vcftools
source activate truvari
ref=/BIGDATA2/scau_jyang_1/pangenome/reference/W64_addSscrofa11.1ChrY.primary_assembly.fa
input=/BIGDATA2/scau_jyang_1/pangenome/projects/2023_PigGraphgenome_phase1/2.graph_genome_construct/scratch11
output=/BIGDATA2/scau_jyang_1/pangenome/projects/2023_PigGraphgenome_phase1/2.graph_genome_construct/scratch11

bcftools concat --threads 24 --allow-overlaps \
${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_merge_0.8.vcf.gz \
${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_merge_0.8.vcf.gz \
${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_merge_0.8.vcf.gz \
${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_merge_0.8.vcf.gz | \
bcftools sort --temp-dir ${output} -O z -o ${output}/6.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8.vcf.gz
tabix -f -p vcf ${output}/6.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8.vcf.gz

