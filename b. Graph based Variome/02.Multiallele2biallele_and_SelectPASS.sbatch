#!/bin/bash
# usage:sbatch -a 0-96 2.Multiallele2biallele_and_SelectPASS.sbatch 1 97
# Job name:
#SBATCH --job-name=2.Multiallele2biallele_and_SelectPASS
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=24
# Memory per node:
#SBATCH --mem=58G
#SBATCH --output=2.Multiallele2biallele_and_SelectPASS_%A_%a_%j.out
start=$1
end=$2
input=./
output=./
mkdir -p ${output}
array=(`cat ${input}/5.totalsample_arry_file.txt | awk '{if(NR>='"$start"'&&NR<='"$end"')print}'`)

module load anaconda3
module load parallel
source activate samtools
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C

echo $(date +"%Y-%m-%d %T") begin split and PASS ${array[SLURM_ARRAY_TASK_ID]}
parallel -j5 -a ${input}/5.totalsample_array/${array[SLURM_ARRAY_TASK_ID]} --memfree 8G " bcftools norm -m-any --multi-overlaps . ${output}/{}.vcf.gz | bcftools view -i 'FILTER=\"PASS\"' -Oz -o ${output}/{}.split.PASS.vcf.gz"
echo $(date +"%Y-%m-%d %T") done split ${array[SLURM_ARRAY_TASK_ID]}

echo $(date +"%Y-%m-%d %T") begin tabix ${array[SLURM_ARRAY_TASK_ID]}
parallel -j5 -a ${input}/5.totalsample_array/${array[SLURM_ARRAY_TASK_ID]} --memfree 8G "tabix -p vcf -f ${output}/{}.split.PASS.vcf.gz"
echo $(date +"%Y-%m-%d %T") done tabix ${array[SLURM_ARRAY_TASK_ID]}