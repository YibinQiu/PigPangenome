#!/bin/bash
# usage:sbatch 4.PanSV_QC.sbatch
# Job name:
#SBATCH --job-name=4.PanSV_QC
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-low-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=100MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/4.PanSV_QC_%j.out
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

plink --vcf ${input}/PigPanSV.vcf.gz --allow-extra-chr --make-bed --out ${output}/PigPanSV_autosome --chr 1-18 --keep-allele-order --threads 64 --vcf-half-call missing --geno 0.1
echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID