#!/bin/bash
# usage:sbatch --export=pop1="DUR",pop2="20MDUR" 11.Fst_WENS_vs_USA.sbatch
# usage:sbatch --export=pop1="LR",pop2="20MLR" 11.Fst_WENS_vs_USA.sbatch
# usage:sbatch --export=pop1="LW",pop2="20MLW" 11.Fst_WENS_vs_USA.sbatch
# Job name:
#SBATCH --job-name=11.Fst_WENS_vs_USA
#SBATCH --account=szzx
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000MB
#SBATCH -t 48:00:00
#SBATCH --array=0-17%18
#SBATCH --output=logs/11.Fst_WENS_vs_USA_%A_%a_%j.out
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
source ~/software/XP-CLR_Env/bin/activate
input=./
output=./
group=${input}/group_1
chr=(`seq 1 18`)

if [ ! -d ${output}/${pop1}_vs_${pop2}_SNP ]
then mkdir -p ${output}/${pop1}_vs_${pop2}_SNP
fi

if [ ! -d ${output}/${pop1}_vs_${pop2}_Indel ]
then mkdir -p ${output}/${pop1}_vs_${pop2}_Indel
fi

if [ ! -d ${output}/${pop1}_vs_${pop2}_SV ]
then mkdir -p ${output}/${pop1}_vs_${pop2}_SV
fi

echo  $(date +"%Y-%m-%d %T") start ${pop1}_${pop2}

vcftools --gzvcf ${input}/PigPanSNP_autosome.vcf.gz --weir-fst-pop ${group}/${pop1}.pop --weir-fst-pop ${group}/${pop2}.pop --chr chr${chr[SLURM_ARRAY_TASK_ID]}  --fst-window-size 50000 --fst-window-step 25000 --out ${output}/${pop1}_vs_${pop2}_SNP/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst > ${output}/${pop1}_vs_${pop2}_SNP/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst.log 2>&1

vcftools --gzvcf ${input}/PigPanIndel_autosome.vcf.gz --weir-fst-pop ${group}/${pop1}.pop --weir-fst-pop ${group}/${pop2}.pop --chr chr${chr[SLURM_ARRAY_TASK_ID]}  --fst-window-size 50000 --fst-window-step 25000 --out ${output}/${pop1}_vs_${pop2}_Indel/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst > ${output}/${pop1}_vs_${pop2}_Indel/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst.log 2>&1

vcftools --gzvcf ${input}/PigPanSV_autosome.vcf.gz --weir-fst-pop ${group}/${pop1}.pop --weir-fst-pop ${group}/${pop2}.pop --chr chr${chr[SLURM_ARRAY_TASK_ID]}  --fst-window-size 50000 --fst-window-step 25000 --out ${output}/${pop1}_vs_${pop2}_SV/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst > ${output}/${pop1}_vs_${pop2}_SV/${pop1}_${pop2}_chr${chr[SLURM_ARRAY_TASK_ID]}.fst.log 2>&1

echo  $(date +"%Y-%m-%d %T") done ${pop1}_${pop2} 

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID