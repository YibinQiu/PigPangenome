#!/bin/bash
# usage:sbatch 07.QTLAnnotation_Fst.sbatch
# Job name:
#SBATCH --job-name=07.QTLAnnotation_Fst
#SBATCH --account=szzx
#SBATCH --qos=nolimit
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=gpu
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1500MB
#SBATCH --output=logs/07.QTLAnnotation_Fst_%j.out
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
mkdir -p ~/RTEMP
export TMPDIR=~/RTEMP
script=./

echo  $(date +"%Y-%m-%d %T") start
input=./DUR_Gen
cd ${input}
Rscript ${script}/QTLAnnotation_Fst.R DUR ADUR 5000 > variCombinedQTLAnnotation_R.log 2>&1

input=./LR_Gen
cd ${input}
Rscript ${script}/QTLAnnotation_Fst.R LR WLR 5000 > variCombinedQTLAnnotation_R.log 2>&1
echo  $(date +"%Y-%m-%d %T") all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID