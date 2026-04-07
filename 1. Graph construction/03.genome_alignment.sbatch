#!/bin/bash
# usage:sbatch 3.genome_alignment.sbatch
#SBATCH --job-name=3.genome_alignment
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --qos=nolimit
#SBATCH --partition=cpu-high-mem,fat,fat-8253
#SBATCH --cpus-per-task=20
#SBATCH --mem=930G
#SBATCH --array=0-31%16
#SBATCH --output=logs/3.genome_alignment_%A_%a_%j.out
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
mkdir -p ${output}/${simple[SLURM_ARRAY_TASK_ID]}

echo $(date +"%m-%d-%Y-%T"): start whole genome alignment for ${genome[SLURM_ARRAY_TASK_ID]}:${simple[SLURM_ARRAY_TASK_ID]} alignment

#parallel -a ${output}/${simple[SLURM_ARRAY_TASK_ID]}/seqID.txt -j 4 --joblog ${output}/${simple[SLURM_ARRAY_TASK_ID]}/jobs.log "/usr/bin/time -v nucmer --maxmatch -c 100 -b 500 -l 50 -t 2 -p ${output}/${simple[SLURM_ARRAY_TASK_ID]}/{} ${ref} ${output}/${simple[SLURM_ARRAY_TASK_ID]}/{}.fa"      # Whole genome alignment.

nucmer --mum -c 1000 --maxgap=500 -t 20 -p ${output}/${simple[SLURM_ARRAY_TASK_ID]} ${ref}/W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa ${input}/${genome[SLURM_ARRAY_TASK_ID]}
delta-filter -i 90 -l 1000 -1 ${output}/${simple[SLURM_ARRAY_TASK_ID]}.delta > ${output}/${simple[SLURM_ARRAY_TASK_ID]}.filter.delta

echo $(date +"%m-%d-%Y-%T"): done for ${genome[SLURM_ARRAY_TASK_ID]} alignment
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID