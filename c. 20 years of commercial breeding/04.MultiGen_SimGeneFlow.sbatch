#!/bin/bash
# Job name:
# usage:sbatch 04.MultiGen_SimGeneFlow.sbatch
#SBATCH --job-name=04.MultiGen_SimGeneFlow
#SBATCH --qos=nolimit
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=1300MB
#SBATCH --array 0-999%30
#SBATCH --output=04.MultiGen_SimGeneFlow_%A_%a_%j.out
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
#for i in `seq 1 1000`; do cp MultiGen_SimGeneFlow_test.prm MultiGen_SimGeneFlow_test${i}.prm; sed -i "s/Multi Generation Selection Simulation/Multi Generation Selection Simulation ${i}/g" MultiGen_SimGeneFlow_test${i}.prm; done
input=./

rep=(`seq 1 1000`)
echo $(date +"%Y-%m-%d %T") QMsim for MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}S21.prm
echo "++++++++++++++++++++++++++++++++++++++++"
mkdir -p ${input}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}_S21
output=${input}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}_S21
cp ${input}/QMSim ${output}/
cp ${input}/MultiGen_SimGeneFlow_S21.prm ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm
sed -i "s/Multi Generation Selection Simulation/Multi Generation Selection Simulation ${rep[SLURM_ARRAY_TASK_ID]}/g" ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm

sed -i "s#qiuyb#${input}/r_MultiGen_SimGeneFlow_test${rep[SLURM_ARRAY_TASK_ID]}/seed#g" ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm

cd ${output}
./QMSim ./MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm -o -16bit
#######################################################
#######################################################
input=./
rep=(`seq 1 1000`)
echo $(date +"%Y-%m-%d %T") QMsim for MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}W51.prm
echo "++++++++++++++++++++++++++++++++++++++++"
mkdir -p ${input}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}_W51
output=${input}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}_W51
cp ${input}/QMSim ${output}/
cp ${input}/MultiGen_SimGeneFlow_W51.prm ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm
sed -i "s/Multi Generation Selection Simulation/Multi Generation Selection Simulation ${rep[SLURM_ARRAY_TASK_ID]}/g" ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm

sed -i "s#qiuyb#${input}/r_MultiGen_SimGeneFlow2_test${rep[SLURM_ARRAY_TASK_ID]}/seed#g" ${output}/MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm

cd ${output}
./QMSim ./MultiGen_SimGeneFlow_${rep[SLURM_ARRAY_TASK_ID]}.prm -o -16bit
echo "++++++++++++++++++++++++++++++++++++++++"
echo $(date +"%Y-%m-%d %T") done for MultiGen_SimGeneFlow_test${rep[SLURM_ARRAY_TASK_ID]}.prm
echo "processs will sleep 10s"
sleep 10
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID