#!/bin/bash
# usage:sbatch --export=pop="DUR",pop1="03ADUR",pop2="13ADUR" 01.PigPan_Fst.sbatch
# usage:sbatch --export=pop="DUR",pop1="03ADUR",pop2="21ADUR" 01.PigPan_Fst.sbatch
# usage:sbatch --export=pop="DUR",pop1="13ADUR",pop2="21ADUR" 01.PigPan_Fst.sbatch

# usage:sbatch --export=pop="LR",pop1="03WLR",pop2="12WLR" 01.PigPan_Fst.sbatch
# usage:sbatch --export=pop="LR",pop1="03WLR",pop2="18WLR" 01.PigPan_Fst.sbatch
# usage:sbatch --export=pop="LR",pop1="12WLR",pop2="18WLR" 01.PigPan_Fst.sbatch
# Job name:
#SBATCH --job-name=01.PigPan_Fst
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=15000MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/01.PigPan_Fst_%j.out
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
conda activate tools_env

input=./
output=./

mkdir -p ${output}/${pop1}_${pop2}
gcta --bfile ${input}/PigPanSNP_autosome --fst --geno 0.1 --sub-popu ${input}/${pop1}_${pop2}.subpop --out ${output}/${pop1}_${pop2}_SNP --threads 64
gcta --bfile ${input}/PigPanIndel_autosome --fst --geno 0.1 --sub-popu ${input}/${pop1}_${pop2}.subpop --out ${output}/${pop1}_${pop2}_Indel --threads 64
gcta --bfile ${input}/PigPanSV_autosome --fst --geno 0.1 --sub-popu ${input}/${pop1}_${pop2}.subpop --out ${output}/${pop1}_${pop2}_SV --threads 64
echo  $(date +"%Y-%m-%d %T")  done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID