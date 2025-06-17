#!/bin/bash
# usage:sbatch 03.Caculate_Fst_Ref_freq.sbatch
# Job name:
#SBATCH --job-name=03.Caculate_Fst_Ref_freq
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15000MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/03.Caculate_Fst_Ref_freq_%j.out
CURDIR=`pwd`
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID
NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`
for i in $NODES
do
echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID
done

echo "process will start at : "
date
echo "++++++++++++++++++++++++++++++++++++++++"
. /public/home/szzx/qiuyb/miniconda3/bin/activate
conda activate tools_env
input2=./
input1=./
output=./

mkdir -p ${output}/DUR_Gen

##for SNP 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
rm ${output}/DUR_Gen/fst_snp01.ID ${output}/DUR_Gen/fst_Indel01.ID ${output}/DUR_Gen/fst_sv01.ID
echo $(date +"%Y-%m-%d %T") for SNP 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
cat ${output}/4.PigPanSNP_Fst/03ADUR/03ADUR_13ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_snp01.ID
#cat ${output}/4.PigPanSNP_Fst/03ADUR/03ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_snp01.ID
cat ${output}/4.PigPanSNP_Fst/13ADUR/13ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_snp01.ID

##for Indel 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
echo $(date +"%Y-%m-%d %T") for Indel 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
cat ${output}/4.PigPanIndel_Fst/03ADUR/03ADUR_13ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_Indel01.ID
#cat ${output}/4.PigPanIndel_Fst/03ADUR/03ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_Indel01.ID
cat ${output}/4.PigPanIndel_Fst/13ADUR/13ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_Indel01.ID

##for SV 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
echo $(date +"%Y-%m-%d %T") for SV 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
cat ${output}/3.PigPanSV_Fst/03ADUR/03ADUR_13ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_sv01.ID
#cat ${output}/3.PigPanSV_Fst/03ADUR/03ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_sv01.ID
cat ${output}/3.PigPanSV_Fst/13ADUR/13ADUR_21ADUR.fst_0.1.ID >> ${output}/DUR_Gen/fst_sv01.ID

###################################################################################################
echo $(date +"%Y-%m-%d %T") make cluster and keepSample.txt
rm ${output}/DUR_Gen/cluster.txt
cat ${output}/cluster.txt |grep -P "03ADUR|13ADUR|21ADUR" >>${output}/DUR_Gen/cluster.txt
cat ${output}/DUR_Gen/cluster.txt |awk '{print $1}' > ${output}/DUR_Gen/keepSample.txt

echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst SNP position and refmaf
vcftools --gzvcf ${input1}/4.PanSNP_autosome.vcf.gz --keep ${output}/DUR_Gen/keepSample.txt --snps ${output}/DUR_Gen/fst_snp01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/DUR_Gen/PigPanSNP_ADUR.vcf.gz

plink --vcf ${output}/DUR_Gen/PigPanSNP_ADUR.vcf.gz --keep-allele-order --make-bed --out ${output}/DUR_Gen/PigPanSNP_ADUR --vcf-half-call missing
plink --bfile ${output}/DUR_Gen/PigPanSNP_ADUR --keep-allele-order --freq --within ${output}/DUR_Gen/cluster.txt --out ${output}/DUR_Gen/ADUR_snp01freq.txt
###
echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst Indel IDition and refmaf
vcftools --gzvcf ${input1}/4.PanIndel_autosome.vcf.gz --keep ${output}/DUR_Gen/keepSample.txt --snps ${output}/DUR_Gen/fst_Indel01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/DUR_Gen/PigPanIndel_ADUR.vcf.gz

plink --vcf ${output}/DUR_Gen/PigPanIndel_ADUR.vcf.gz --keep-allele-order --make-bed --out ${output}/DUR_Gen/PigPanIndel_ADUR --vcf-half-call missing
plink --bfile ${output}/DUR_Gen/PigPanIndel_ADUR --keep-allele-order --freq --within ${output}/DUR_Gen/cluster.txt --out ${output}/DUR_Gen/ADUR_Indel01freq.txt
###
echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst SV IDition and refmaf
vcftools --gzvcf ${input2}/PanSV_autosome.vcf.gz --keep ${output}/DUR_Gen/keepSample.txt --snps ${output}/DUR_Gen/fst_sv01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/DUR_Gen/PigPanSV_ADUR.vcf.gz

plink --vcf ${output}/DUR_Gen/PigPanSV_ADUR.vcf.gz --keep-allele-order --make-bed --out ${output}/DUR_Gen/PigPanSV_ADUR --vcf-half-call missing
plink --bfile ${output}/DUR_Gen/PigPanSV_ADUR --keep-allele-order --freq --within ${output}/DUR_Gen/cluster.txt --out ${output}/DUR_Gen/ADUR_sv01freq.txt

##################################################################################################
mkdir -p ${output}/LR_Gen

##for SNP 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
rm ${output}/LR_Gen/fst_snp01.ID ${output}/LR_Gen/fst_indel01.ID ${output}/LR_Gen/fst_sv01.ID
echo $(date +"%Y-%m-%d %T") for SNP 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR

cat ${output}/4.PigPanSNP_Fst/03WLR/03WLR_12WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_snp01.ID
#cat ${output}/4.PigPanSNP_Fst/03WLR/03WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_snp01.ID
cat ${output}/4.PigPanSNP_Fst/12WLR/12WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_snp01.ID

##for INDEL 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
echo $(date +"%Y-%m-%d %T") for INDEL 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR

cat ${output}/4.PigPanINDEL_Fst/03WLR/03WLR_12WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_indel01.ID
#cat ${output}/4.PigPanINDEL_Fst/03WLR/03WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_indel01.ID
cat ${output}/4.PigPanINDEL_Fst/12WLR/12WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_indel01.ID

##for SV 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR
echo $(date +"%Y-%m-%d %T") for SV 03ADUR 13ADUR 21ADUR 03WLR 12WLR 18WLR

cat ${output}/3.PigPanSV_Fst/03WLR/03WLR_12WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_sv01.ID
#cat ${output}/3.PigPanSV_Fst/03WLR/03WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_sv01.ID
cat ${output}/3.PigPanSV_Fst/12WLR/12WLR_18WLR.fst_0.1.ID >> ${output}/LR_Gen/fst_sv01.ID

##################################################################################################
echo $(date +"%Y-%m-%d %T") make cluster and keepSample.txt
rm ${output}/LR_Gen/cluster.txt
cat ${output}/cluster.txt |grep -P "03WLR|12WLR|18WLR" >>${output}/LR_Gen/cluster.txt
cat ${output}/LR_Gen/cluster.txt |awk '{print $1}' > ${output}/LR_Gen/keepSample.txt
echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst SNP position and refmaf

vcftools --gzvcf ${input1}/4.PanSNP_autosome.vcf.gz --keep ${output}/LR_Gen/keepSample.txt --snps ${output}/LR_Gen/fst_snp01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/LR_Gen/PigPanSNP_WLR.vcf.gz

plink --vcf ${output}/LR_Gen/PigPanSNP_WLR.vcf.gz --keep-allele-order --make-bed --out ${output}/LR_Gen/PigPanSNP_WLR --vcf-half-call missing
plink --bfile ${output}/LR_Gen/PigPanSNP_WLR --keep-allele-order --freq --within ${output}/LR_Gen/cluster.txt --out ${output}/LR_Gen/WLR_snp01freq.txt
###
echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst Indel IDition and refmaf
vcftools --gzvcf ${input1}/4.PanIndel_autosome.vcf.gz --keep ${output}/LR_Gen/keepSample.txt --snps ${output}/LR_Gen/fst_Indel01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/LR_Gen/PigPanIndel_WLR.vcf.gz

plink --vcf ${output}/LR_Gen/PigPanIndel_WLR.vcf.gz --keep-allele-order --make-bed --out ${output}/LR_Gen/PigPanIndel_WLR --vcf-half-call missing
plink --bfile ${output}/LR_Gen/PigPanIndel_WLR --keep-allele-order --freq --within ${output}/LR_Gen/cluster.txt --out ${output}/LR_Gen/WLR_Indel01freq.txt
###
echo $(date +"%Y-%m-%d %T") keep keepSample.txt and fst SV IDition and refmaf
vcftools --gzvcf ${input2}/PanSV_autosome.vcf.gz --keep ${output}/LR_Gen/keepSample.txt --snps ${output}/LR_Gen/fst_sv01.ID --recode --recode-INFO-all --stdout | bgzip -c > ${output}/LR_Gen/PigPanSV_WLR.vcf.gz

plink --vcf ${output}/LR_Gen/PigPanSV_WLR.vcf.gz --keep-allele-order --make-bed --out ${output}/LR_Gen/PigPanSV_WLR --vcf-half-call missing
plink --bfile ${output}/LR_Gen/PigPanSV_WLR --keep-allele-order --freq --within ${output}/LR_Gen/cluster.txt --out ${output}/LR_Gen/WLR_sv01freq.txt
echo  $(date +"%Y-%m-%d %T") all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID