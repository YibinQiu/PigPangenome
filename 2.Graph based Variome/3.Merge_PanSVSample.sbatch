#!/bin/bash
# usage:sbatch 3.Merge_PanSVSample.sbatch
# Job name:
#SBATCH --job-name=3.Merge_PanSVSample.sbatch
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
#SBATCH --output=3.Merge_PanSVSample_%j.out

. /public/home/szzx/qiuyb/miniconda3/bin/activate
conda activate truvari_env

ls ${output}/7.totalsample_batch_file/totalsample_batch_*.vcf.gz > ${output}/9.merge_totalsample_array.list
echo $(date +"%Y-%m-%d %T") begin merge_batch
bcftools merge -m none -l ${output}/9.merge_totalsample_array.list --threads 23 -Oz -o ${output}/9.merge_totalsample_array.vcf.gz
tabix -p vcf -f ${output}/9.merge_totalsample_array.vcf.gz


bcftools sort -T ${output} ${output}/9.merge_totalsample_array.vcf.gz -O z -o ${output}/9.merge_totalsample_array_norm_sort.vcf.gz

echo $(date +"%Y-%m-%d %T") begin anno
zcat ${output}/9.merge_totalsample_array_norm_sort.vcf.gz | truvari anno svinfo -m 10 | bgzip -f --threads 64 > ${output}/9.merge_totalsample_array_norm_sort_anno.vcf.gz
tabix -f -p vcf ${output}/9.merge_totalsample_array_norm_sort_anno.vcf.gz

echo $(date +"%Y-%m-%d %T") begin keep DEL and UNK
bcftools filter -i 'SVTYPE="DEL" || SVTYPE="UNK"' ${output}/9.merge_totalsample_array_norm_sort_anno.vcf.gz -Oz -o ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK.vcf.gz
tabix -f -p vcf ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK.vcf.gz

echo $(date +"%Y-%m-%d %T") begin keep INS
bcftools filter -i 'SVTYPE="INS"' ${output}/9.merge_totalsample_array_norm_sort_anno.vcf.gz -Oz -o ${output}/9.merge_totalsample_array_norm_sort_anno_INS.vcf.gz
tabix -f -p vcf ${output}/9.merge_totalsample_array_norm_sort_anno_INS.vcf.gz

echo $(date +"%Y-%m-%d %T") begin merge DEL and UNK
truvari collapse -i ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK.vcf.gz -o ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_merge_0.8.vcf -c ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_collapse_0.8.vcf -f $ref -k common -r 1000 -p 0 -B 30 -P 0.0 -O 0.8 -s 20 -S 5000000
bgzip --threads 24 -f ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_merge_0.8.vcf.gz
bgzip --threads 24 -f ${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_collapse_0.8.vcf

echo $(date +"%Y-%m-%d %T") begin merge INS
truvari collapse -i ${output}/9.merge_totalsample_array_norm_sort_anno_INS.vcf.gz -o ${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_merge_0.8.vcf -c ${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_collapse_0.8.vcf -f $ref -k common -r 20 -p 0 -B 30 -P 0.8 -s 20 -S 5000000
bgzip --threads 24 -f ${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_merge_0.8.vcf
tabix -f -p vcf ${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_merge_0.8.vcf.gz
bgzip --threads 24 -f ${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_collapse_0.8.vcf

echo $(date +"%Y-%m-%d %T") begin merge all
bcftools concat --threads 64 --allow-overlaps \
${output}/9.merge_totalsample_array_norm_sort_anno_DEL_UNK_truvari_merge_0.8.vcf.gz \
${output}/9.merge_totalsample_array_norm_sort_anno_INS_truvari_merge_0.8.vcf.gz | \
bcftools sort --temp-dir ${output} -O z -o ${output}/PanSV.vcf.gz
tabix -f -p vcf ${output}/PanSV.vcf.gz