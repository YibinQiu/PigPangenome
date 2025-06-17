#!/bin/bash
# usage:sbatch --export=epi_array_task="epi_array_task" --array=0-60%61 02.epistasis.sbatch

#SBATCH --job-name=02.epistasis
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --qos=nolimit
#SBATCH --partition=cpu-high-mem
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=15000MB
#SBATCH --output=logs/02.epistasis_%A_%a_%j.out
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

task_trait=(`cat ${output}/${epi_array_task} |awk '{print $1}'`)
task_set=(`cat ${output}/${epi_array_task} |awk '{print $2}'`)
task_phe=(`cat ${output}/${epi_array_task} |awk '{print $3}'`)

echo $(date +"%Y-%m-%d %T") start for ${task_trait[SLURM_ARRAY_TASK_ID]} ${task_set[SLURM_ARRAY_TASK_ID]} ${task_phe[SLURM_ARRAY_TASK_ID]}
#mkdir -p ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}
plink --allow-no-sex \
       --bfile ${input}/GWASsample_withDUP \
       --pheno ${input}/adjusted_phe/${task_phe[SLURM_ARRAY_TASK_ID]} \
       --epistasis set-by-all \
       --set ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/${task_set[SLURM_ARRAY_TASK_ID]} \
       --epi1 0.000001 \
       --threads 30\
       --noweb \
       --out ${output}/${task_trait[SLURM_ARRAY_TASK_ID]}/${task_trait[SLURM_ARRAY_TASK_ID]}
     
echo $(date +"%m-%d-%Y-%T"): done for ${task_trait[SLURM_ARRAY_TASK_ID]} ${task_set[SLURM_ARRAY_TASK_ID]} ${task_phe[SLURM_ARRAY_TASK_ID]}
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID