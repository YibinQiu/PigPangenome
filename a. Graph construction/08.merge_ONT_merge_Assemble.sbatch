#!/bin/bash
# usage:yhbatch 8.merge_ONT_merge_Assemble.sbatch
# Job name:
#SBATCH --job-name=8.merge_ONT_merge_Assemble.sbatch
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
#SBATCH --output=8.merge_ONT_merge_Assemble_%j.out

module load bcftools

input=/BIGDATA2/scau_jyang_1/pangenome/projects/2023_PigGraphgenome_phase1/2.graph_genome_construct/scratch11
output=/BIGDATA2/scau_jyang_1/pangenome/projects/2023_PigGraphgenome_phase1/2.graph_genome_construct/scratch11

ontsamples=(`cat ${output}/0.ONT_sample`)

for ontsample in ${ontsamples[*]}
do
    echo ${input}/0.ONT_highConfidence/${ontsample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz >> ${output}/2.ONT_sample_path.list
done
echo $(date +"%Y-%m-%d %T") start ONT_merge
bcftools merge -m none --file-list ${output}/2.ONT_sample_path.list -O z -o ${output}/2.ONT_merge.vcf.gz
tabix -f -p vcf ${output}/2.ONT_merge.vcf.gz

Chromosomesamples=(`cat ${output}/0.Chromosome_sample`)
for Chromosomesample in ${Chromosomesamples[*]}
do
    echo ${input}/0.Assemble_highConfidence/${Chromosomesample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz >> ${output}/2.Chromosome_Scaffold_Contig_sample_path.list
done

Scaffold_Contigsamples=(`cat ${output}/0.Scaffold_Contig_sample`)
for Scaffold_Contigsample in ${Scaffold_Contigsamples[*]}
do
    echo ${input}/0.Assemble_highConfidence/${Scaffold_Contigsample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz >> ${output}/2.Chromosome_Scaffold_Contig_sample_path.list
done
echo $(date +"%Y-%m-%d %T") start Chromosome_Scaffold_Contig_merge
bcftools merge -m none --file-list ${output}/2.Chromosome_Scaffold_Contig_sample_path.list -O z -o ${output}/2.Chromosome_Scaffold_Contig_merge.vcf.gz
tabix -f -p vcf ${output}/2.Chromosome_Scaffold_Contig_merge.vcf.gz
echo $(date +"%Y-%m-%d %T") all done

