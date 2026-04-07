#!/bin/bash
# usage:sbatch --export=epi_array_task="epi_array_task" --array=0-60%61 05.pick_tag_epi.sbatch
# Job name:
#SBATCH --job-name=05.pick_tag_epi
#SBATCH --partition=fat-8253,cpu-low-mem,cpu-high-mem,cpu-6226
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --account=szzx
# Processors per task:
#SBATCH --cpus-per-task=64
#
# Memory per node:
#SBATCH --mem-per-cpu=7500MB
#SBATCH --output=logs/05.pick_tag_epi_%A_%a_%j.out

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
conda activate tools_env


epi=./
chr=./
output=./
task_trait=(`cat ${output}/${epi_array_task} |awk '{print $1}'`)

echo $(date +"%Y-%m-%d %T") start for ${task_trait[SLURM_ARRAY_TASK_ID]}
mkdir -p ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/multithreads
python3 7.pick_epi.py --trait ${task_trait[SLURM_ARRAY_TASK_ID]} --pvalue "1e-8" --rsqure 0.1 --chrLen ${chr} --epi ${epi}/${task_trait[SLURM_ARRAY_TASK_ID]} --output ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/multithreads

echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 60
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID