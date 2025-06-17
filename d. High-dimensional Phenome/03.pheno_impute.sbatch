#!/bin/bash
#usage: sbatch 03.pheno_impute.sbatch
#SBATCH --job-name=03.pheno_impute
#SBATCH --qos=nolimit
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --output=logs/03.pheno_impute_%A_%a.out
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=7500MB
#SBATCH --cpus-per-task=5
#SBATCH --partition=cpu-high-mem,cpu-low-mem,cpu-6226
#SBATCH --time=148:00:00

# ============================================================
CURDIR=`pwd`
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID
NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`
for i in $NODES
do
echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID
done
# Load required modules
. ~/miniconda3/bin/activate
conda activate ac
# fit.py can be obtained from https://github.com/sriramlab/AutoComplete
python3.7 ~/software/AutoComplete/fit.py \
 datasets/phenotypes/data.csv \
 --id_name ID \
 --copymask_amount 0.5 \
 --batch_size 2048 \
 --epochs 500 --lr 0.1 --quality

echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 60
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID