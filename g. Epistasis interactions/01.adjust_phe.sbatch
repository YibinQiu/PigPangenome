#!/bin/bash
# usage:sbatch 01.adjust_phe.sbatch
# Job name:
#SBATCH --job-name=01.adjust_phe
#SBATCH --partition=cpu-high-mem,cpu-6226,cpu-low-mem
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --account=szzx
# Processors per task:
#SBATCH --cpus-per-task=10
#
# Memory per node:
#SBATCH --mem-per-cpu=7500MB
#SBATCH --output=logs/01.adjust_phe_%j.out

# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=148:00:00

# ============================================================
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
. ~/miniconda3/bin/activate
conda activate r_env

onlysigfactor="FALSE"
input_export="./"
input_pca="./"
input_qcov="./"
input_cov="./"
input_trait="./"

export_list="export_list.txt"
pcafile="S21_S22_S23_merge_GRM_PCA.eigenvec"
output="./"
mkdir -p ${output}
echo $(date +"%m-%d-%Y-%T"): adjust phe start

for trait in `cat ${input_export}/${export_list}| awk '{print $1}'|sort|uniq`
do
  Rscript adjust_phe.R --onlysigfactor ${onlysigfactor} --input_export ${input_export} --input_pca ${input_pca} --input_qcov ${input_qcov} --input_cov ${input_cov} --input_trait ${input_trait} --export_list ${export_list} --pcafile ${pcafile} --traitfile ${trait} --output ${output}
done

echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 60
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID