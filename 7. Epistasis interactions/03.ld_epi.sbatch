#!/bin/bash
# usage:sbatch --export=epi_array_task="epi_array_task" --array=0-60%61 03.ld_epi.sbatch

#SBATCH --job-name=03.ld_epi
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --qos=nolimit
#SBATCH --partition=cpu-high-mem,cpu-low-mem,cpu-6226,fat,fat-8253,gpu
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=7500MB
#SBATCH --output=logs/03.ld_epi_%A_%a_%j.out
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
input=./
output=./
mkdir -p ${output}

task_trait=(`cat ${output}/${epi_array_task} |awk '{print $1}'`)
task_set=(`cat ${output}/${epi_array_task} |awk '{print $2}'`)

plink --bfile ${input}/0.data/GWASsample_withDUP --ld-snp-list ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/${task_set[SLURM_ARRAY_TASK_ID]} --ld-window 999999999 --ld-window-kb 2000 --ld-window-r2 0.1 --r2 'yes-really' --out ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/pval_0.001_ld --threads 32

cat ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/pval_0.001_ld.ld | awk '$3!=$6{print $3,$6,$7}' > ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/pval_0.001_ld.ld2
mv ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/pval_0.001_ld.ld2 ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/pval_0.001_ld.ld
echo $(date +"%Y-%m-%d %T") start for ld


echo $(date +"%m-%d-%Y-%T"): all done 
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID