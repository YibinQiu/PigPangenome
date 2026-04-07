#!/bin/bash
# usage:sbatch --export=pop1="WS",pop2="MJB" 12.QTLAnnotation_WENS_vs_USA.sbatch
# Job name:
#SBATCH --job-name=12.QTLAnnotation_XPCLR
#SBATCH --account=szzx
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=13000MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/12.QTLAnnotation_XPCLR_%j.out
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
conda activate r_env
script=./
output=./

cd ${output}
Rscript ${script}/QTLAnnotation_WENS_vs_USA.R ${pop1} ${pop2} 5000 > ${pop1}_${pop2}_QTLAnnotation_5000bp.log 2>&1 
echo  $(date +"%Y-%m-%d %T") done ${pop1}_${pop2} QTL enrich

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID