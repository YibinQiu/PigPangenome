#!/bin/bash
# usage:sbatch 07.PigPan_PCA.sbatch
# Job name:
#SBATCH --job-name=7.PigPan_PCA
#SBATCH --account=szzx
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=15000MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/7.PigPan_PCA_%j.out
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
input=./
output=./

#####SNP
plink --bfile ${output}/PanSNP_autosome -keep-allele-order --maf 0.05 --geno 0.1 --make-bed --out ${output}/PanSNP_autosome_filter --allow-extra-chr
gcta --make-grm --out ${output}/PanSNP_autosome_filter_GRM --bfile ${output}/PanSNP_autosome_filter --threads 64
gcta --grm ${output}/PanSNP_autosome_filter_GRM --pca 50 --out ${output}/PanSNP_autosome_filter_PCA --threads 64

#####INDEL
plink --bfile ${output}/PanIndel_autosome -keep-allele-order --maf 0.05 --geno 0.1 --make-bed --out ${output}/PanIndel_autosome_filter --allow-extra-chr
gcta --make-grm --out ${output}/PanIndel_autosome_filter_GRM --bfile ${output}/PanIndel_autosome_filter --threads 64
gcta --grm ${output}/PanIndel_autosome_filter_GRM --pca 50 --out ${output}/PanIndel_autosome_filter_PCA --threads 64

#####SV
plink --bfile ${output}/PanSV_autosome -keep-allele-order --maf 0.05 --geno 0.1 --make-bed --out ${output}/PanSV_autosome_filter --allow-extra-chr
gcta --make-grm --out ${output}/PanSV_autosome_filter_GRM --bfile ${output}/PanSV_autosome_filter --threads 64
gcta --grm ${output}/PanSV_autosome_filter_GRM --pca 50 --out ${output}/PanSV_autosome_filter_PCA --threads 64


echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID