#!/bin/bash
# usage:sbatch 5.svmu2vcf.sbatch
#SBATCH --job-name=5.svmu2vcf
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --qos=nolimit
#SBATCH --partition=cpu-high-mem,fat,fat-8253
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7500MB
#SBATCH --array=0-31%16
#SBATCH --output=logs/5.svmu2vcf_%A_%a_%j.out
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
conda activate svmtools_env
ref=./
qury_path=./
sv_path=./
output=./

genome=(`cat fastaSimpleName.txt | awk '{print $1}'`)
simple=(`cat fastaSimpleName.txt | awk '{print $2}'`)

echo $(date +"%m-%d-%Y-%T"): start for ${genome[SLURM_ARRAY_TASK_ID]}:${simple[SLURM_ARRAY_TASK_ID]} svmu 2 vcf

Rscript svmu2vcf.R --sample_name ${simple[SLURM_ARRAY_TASK_ID]} --svmu_file ${sv_path}/sv.${simple[SLURM_ARRAY_TASK_ID]}.txt --output_file ${output}/${simple[SLURM_ARRAY_TASK_ID]}.vcf --min_size 30 --max_size 5000000 --ref_fasta ${ref}/W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa --query_fasta ${qury_path}/${genome[SLURM_ARRAY_TASK_ID]} --cores 1 --logging TRUE

echo $(date +"%m-%d-%Y-%T"): done for ${genome[SLURM_ARRAY_TASK_ID]} svmu 2 vcf
echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID