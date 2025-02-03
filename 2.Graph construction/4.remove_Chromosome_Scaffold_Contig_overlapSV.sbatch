#!/bin/bash
# usage:sbatch 4.remove_Chromosome_Scaffold_Contig_overlapSV.sbatch
# Job name:
#SBATCH --job-name=4.remove_Chromosome_Scaffold_Contig_overlapSV.sbatch
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#SBATCH --partition=rhenv
# Processors per task:
#SBATCH --cpus-per-task=24
# Memory per node:
#SBATCH --mem=125G
#SBATCH --output=4.remove_Chromosome_Scaffold_Contig_overlapSV_%j.out

export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C

module load bcftools
module load vcftools
input="./"
output="./"
mkdir -p ${output}

echo $(date +"%Y-%m-%d %T") keep overlap sv name
zcat ${input}/3.ONT_vs_Chromosome_Scaffold_Contig/tp-comp.vcf.gz | grep -v "#" |awk '{print $3}' |awk 'BEGIN {FS=";"} {for(i=1; i<=NF; i++) print $i}' > ${output}/overlapSV.txt

Chromosome_Contigsamples=(`cat ${input}/1.Chromosome_sample`)

for Chromosome_Contigsample in ${Chromosome_Contigsamples[*]}
do
    echo $(date +"%Y-%m-%d %T") rm ${Chromosome_Contigsample} overlap SV
    vcftools --gzvcf ${input}/1.Assemble_highConfidence/${Chromosome_Contigsample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz --exclude ${output}/overlapSV.txt --recode --recode-INFO-all --stdout | bcftools sort | bgzip -f > ${output}/${Chromosome_Contigsample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz
    tabix -f -p vcf ${output}/${Chromosome_Contigsample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz
done
echo $(date +"%Y-%m-%d %T") all done


Scaffold_Contigsamples=(`cat ${input}/1.Scaffold_Contig_sample`)

for Scaffold_Contigsample in ${Scaffold_Contigsamples[*]}
do
    echo $(date +"%Y-%m-%d %T") rm ${Scaffold_Contigsample} overlap SV
    vcftools --gzvcf ${input}/1.Assemble_highConfidence/${Scaffold_Contigsample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz --exclude ${output}/overlapSV.txt --recode --recode-INFO-all --stdout | bcftools sort | bgzip -f > ${output}/${Scaffold_Contigsample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz
    tabix -f -p vcf ${output}/${Scaffold_Contigsample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz
done
echo $(date +"%Y-%m-%d %T") all done
