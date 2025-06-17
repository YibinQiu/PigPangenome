#!/bin/bash
# usage:sbatch --export=samplelist="ontsample" --array=0-149%50 2.sniffles.sbatch
# Job name:
#SBATCH --job-name=2.sniffles
#SBATCH --qos=nolimit
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem,cpu-6226,cpu-low-mem,fat-8253
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=15000MB
#SBATCH --output=2.sniffles_%A_%a_%j.out

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
input="./"
sample=(`cat ${input}/${samplelist}`)
output="./"
ref="./"

echo $(date +"%Y-%m-%d %T") ${sample[SLURM_ARRAY_TASK_ID]} start sniffles
sniffles --allow-overwrite -i ${input}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.bam -v ${output}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf --snf ${output}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.snf --threads 64 --cluster-merge-bnd 50 --minsvlen 30 --output-rnames --reference ${ref}/W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa --max-del-seq-len 100000 --sample-id ${sample[SLURM_ARRAY_TASK_ID]} &> ${output}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}_sniffles.vcf.log
echo $(date +"%Y-%m-%d %T") ${sample[SLURM_ARRAY_TASK_ID]} sniffles done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID