#!/bin/bash
# usage:yhbatch 11.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV.sbatch
# Job name:
#SBATCH --job-name=11.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV.sbatch
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
#SBATCH --output=11.merge_ONT_Chromosome_Scaffold_Contig_RemoveOverlapSV_%j.out
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C

module load anaconda3
module load bcftools
module load vcftools
source activate truvari
ref=W64_addSscrofa11.1ChrY.primary_assembly.fa
input=./
output=./

ontsamples=(`cat ${output}/1.ONT_sample`)

for ontsample in ${ontsamples[*]}
do
    echo ${input}/1.ONT_highConfidence/${ontsample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz >> ${output}/5.ONT_Chromosome_Scaffold_Contig_sample_path.list
done

Chromosomesamples=(`cat ${output}/1.Chromosome_sample`)

for Chromosomesample in ${Chromosomesamples[*]}
do
    echo ${input}/4.Chromosome_Scaffold_Contig_highConfidence_remove_overlapSV/${Chromosomesample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz >> ${output}/5.ONT_Chromosome_Scaffold_Contig_sample_path.list
done

Scaffold_Contigsamples=(`cat ${output}/1.Scaffold_Contig_sample`)

for Scaffold_Contigsample in ${Scaffold_Contigsamples[*]}
do
    echo ${input}/4.Chromosome_Scaffold_Contig_highConfidence_remove_overlapSV/${Scaffold_Contigsample}.highConfidence_normalizeTypes_removeOverlapSV.ism.sorted.vcf.gz >> ${output}/5.ONT_Chromosome_Scaffold_Contig_sample_path.list
done
bcftools merge -m none --threads 24 --file-list ${output}/5.ONT_Chromosome_Scaffold_Contig_sample_path.list -O z -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz

################################################################################################################
bcftools filter -i '((SVTYPE="DEL" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000) || (SVTYPE="INV" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000))' ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz --threads 24 -O z -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB.vcf.gz
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB.vcf.gz

echo $(date +"%Y-%m-%d %T") merge 5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB
truvari collapse -i ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB.vcf.gz -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_merge_0.8.vcf -c ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_collapse_0.8.vcf -f $ref -k common -r 1000 -p 0 -B 30 -P 0.0 -O 0.8 -s 20 -S 5000000
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_merge_0.8.vcf.gz
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_100KB_truvari_collapse_0.8.vcf

################################################################################################################
bcftools filter -i 'SVTYPE="INS" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000' ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz --threads 25 -O z -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB.vcf.gz
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB.vcf.gz

echo $(date +"%Y-%m-%d %T") merge 5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB
truvari collapse -i ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB.vcf.gz -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_merge_0.8.vcf -c ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_collapse_0.8.vcf -f $ref -k common -r 50 -p 0 -B 30 -P 0.8 -s 20 -S 5000000
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_merge_0.8.vcf.gz
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_100KB_truvari_collapse_0.8.vcf

#################################################################################################################
bcftools filter -i '((SVTYPE="DEL" && ABS(SVLEN)>100000) || (SVTYPE="INV" && ABS(SVLEN)>100000))' ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz --threads 24 -O z -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB.vcf.gz
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB.vcf.gz

echo $(date +"%Y-%m-%d %T") merge 5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB
truvari collapse -i ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB.vcf.gz -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_merge_0.8.vcf -c ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_collapse_0.8.vcf -f $ref -k common -r 1000 -p 0 -B 30 -P 0.0 -O 0.8 -s 20 -S 5000000
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_merge_0.8.vcf.gz
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_DEL_INV_larger100KB_truvari_collapse_0.8.vcf

####################################################################################################################
bcftools filter -i 'SVTYPE="INS" && ABS(SVLEN)>100000' ${output}/5.ONT_Chromosome_Scaffold_Contig_merge.vcf.gz --threads 24 -O z -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB.vcf.gz
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB.vcf.gz

echo $(date +"%Y-%m-%d %T") merge 5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB
truvari collapse -i ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB.vcf.gz -o ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_merge_0.8.vcf -c ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_collapse_0.8.vcf -f $ref -k common -r 50 -p 0 -B 30 -P 0.8 -s 20 -S 5000000
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_merge_0.8.vcf.gz
bgzip -f --threads 24 ${output}/5.ONT_Chromosome_Scaffold_Contig_merge_INS_larger100KB_truvari_collapse_0.8.vcf
echo $(date +"%Y-%m-%d %T") merge all done

