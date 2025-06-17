#!/bin/bash
# usage:sbatch --export=pop="DUR",pop1="03ADUR",pop2="13ADUR" 02.Fst_filter.sbatch
# usage:sbatch --export=pop="DUR",pop1="03ADUR",pop2="21ADUR" 02.Fst_filter.sbatch
# usage:sbatch --export=pop="DUR",pop1="13ADUR",pop2="21ADUR" 02.Fst_filter.sbatch

# usage:sbatch --export=pop="LR",pop1="03WLR",pop2="12WLR" 02.Fst_filter.sbatch
# usage:sbatch --export=pop="LR",pop1="03WLR",pop2="18WLR" 02.Fst_filter.sbatch
# usage:sbatch --export=pop="LR",pop1="12WLR",pop2="18WLR" 02.Fst_filter.sbatch
# Job name:
#SBATCH --job-name=02.Fst_filter
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-low-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=7500MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/02.Fst_filter_%j.out
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
conda activate r_env

input=./
output=./

echo  $(date +"%Y-%m-%d %T") fst ${pop}: ${pop1} vs ${pop2} filter

Rscript Fst_filter.R ./ ${pop1} ${pop2} 0.1

Rscript Fst_filter.R ./ ${pop1} ${pop2} 0.1

Rscript Fst_filter.R ./ ${pop1} ${pop2} 0.1

echo  $(date +"%Y-%m-%d %T") all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID