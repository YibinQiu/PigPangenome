#!/bin/bash
# usage:sbatch --export=samplelist="ONTsample.list" --array=1-150%5 1.map.sbatch
# Job name:
#SBATCH --job-name=1.map
#SBATCH --qos=nolimit
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem,cpu-6226,cpu-low-mem,fat-8253
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=7500MB
#SBATCH --output=1.map_%A_%a_%j.out

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
conda activate winnowmap_env
ref="./"
sample=(`cat ${input}/${samplelist}`)
input2="./"
output="./"
temp="./${sample[SLURM_ARRAY_TASK_ID]}"
mkdir -p ${temp}
cd ${output}

echo $(date +"%Y-%m-%d %T") begin ${sample[SLURM_ARRAY_TASK_ID]} W64
winnowmap --MD -t 23 -W ${ref}/W64_repetitive_k15.txt -ax map-ont ${ref}/W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa ${input2}/${sample[SLURM_ARRAY_TASK_ID]}.pass.fastq.gz | samtools sort -@ 23 -m 7G -o ${temp}/${sample[SLURM_ARRAY_TASK_ID]}_W64.bam --write-index -
echo $(date +"%Y-%m-%d %T") ${sample[SLURM_ARRAY_TASK_ID]} done


echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID