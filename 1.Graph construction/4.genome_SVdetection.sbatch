#!/bin/bash
# usage:sbatch 4.genome_SVdetection.sbatch
#SBATCH --job-name=4.genome_SVdetection
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --qos=nolimit
#SBATCH --partition=cpu-high-mem,fat,fat-8253
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6400MB
#SBATCH --array=0-31%16
#SBATCH --output=logs/4.genome_SVdetection_%A_%a_%j.out
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
conda activate svmu_env
ref=./
input=./
output=./

genome=(`cat fastaSimpleName.txt | awk '{print $1}'`)
simple=(`cat fastaSimpleName.txt | awk '{print $2}'`)
#mkdir -p ${output}/${simple[SLURM_ARRAY_TASK_ID]}
cd ${output}
echo $(date +"%m-%d-%Y-%T"): start for ${genome[SLURM_ARRAY_TASK_ID]}:${simple[SLURM_ARRAY_TASK_ID]} svmu

/public/home/szzx/qiuyb/software/svmu/svmu ${output}/${simple[SLURM_ARRAY_TASK_ID]}.filter.delta ${ref}/W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa ${input}/${genome[SLURM_ARRAY_TASK_ID]} l blank.txt ${simple[SLURM_ARRAY_TASK_ID]}

echo $(date +"%m-%d-%Y-%T"): done for ${genome[SLURM_ARRAY_TASK_ID]}:${simple[SLURM_ARRAY_TASK_ID]} svmu
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID