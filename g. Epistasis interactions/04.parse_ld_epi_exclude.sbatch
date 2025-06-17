#!/bin/bash
# usage:sbatch --export=epi_array_task="epi_array_task" --array=0-60%61 04.parse_ld_epi_exclude.sbatch
# Job name:
#SBATCH --job-name=04.parse_ld_epi_exclude
#SBATCH --partition=cpu-high-mem
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --account=szzx
# Processors per task:
#SBATCH --cpus-per-task=24
#
# Memory per node:
#SBATCH --mem-per-cpu=15000MB
#SBATCH --output=logs/04.parse_ld_epi_exclude_%A_%a_%j.out

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

gwas=./
epi=./
ld=./
output=./
include=./merge_SVtransMissing_retain_ID.txt
task_trait=(`cat ${output}/${epi_array_task} |awk '{print $1}'`)

echo $(date +"%Y-%m-%d %T") start for ${task_trait[SLURM_ARRAY_TASK_ID]}

Rscript parse_ld_epi_exclude --trait ${task_trait[SLURM_ARRAY_TASK_ID]} --gwas ${gwas} --epi ${epi}/${task_trait[SLURM_ARRAY_TASK_ID]} --ld ${ld}/${task_trait[SLURM_ARRAY_TASK_ID]} --include ${include} --output ${output}/${task_trait[SLURM_ARRAY_TASK_ID]} 

echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 60
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID