#!/bin/bash
# usage:sbatch --export=input="./DUR_Gen",fst="ADUR_sv01freq.txt.frq.strat",delta1="0.47",delta2="0.59",id1="./PigPanSV_Fst/03ADUR_13ADUR.fst_0.1.ID",id2="./PigPanSV_Fst/13ADUR_21ADUR.fst_0.1.ID",outpre="SV"  06.SelectionPrefernece_Fst.sbatch

# usage:sbatch --export=input="./DUR_Gen",fst="ADUR_snp01freq.txt.frq.strat",delta1="0.47",delta2="0.59",id1="./PigPanSNP_Fst/03ADUR_13ADUR.fst_0.1.ID",id2="./PigPanSNP_Fst/13ADUR_21ADUR.fst_0.1.ID",outpre="SNP"  06.SelectionPrefernece_Fst.sbatch

# usage:sbatch --export=input="./DUR_Gen",fst="ADUR_indel01freq.txt.frq.strat",delta1="0.47",delta2="0.59",id1="./PigPanINDEL_Fst/03ADUR_13ADUR.fst_0.1.ID",id2="./PigPanIndel_Fst/13ADUR_21ADUR.fst_0.1.ID",outpre="INDEL"  06.SelectionPrefernece_Fst.sbatch

# usage:sbatch --export=input="./LR_Gen",fst="WLR_sv01freq.txt.frq.strat",delta1="0.46",delta2="0.44",id1="./PigPanSV_Fst/03WLR_12WLR.fst_0.1.ID",id2="./PigPanSV_Fst/12WLR_18WLR.fst_0.1.ID",outpre="SV"  06.SelectionPrefernece_Fst.sbatch

# usage:sbatch --export=input="./LR_Gen",fst="WLR_snp01freq.txt.frq.strat",delta1="0.46",delta2="0.44",id1="./PigPanSNP_Fst/03WLR_12WLR.fst_0.1.ID",id2="./PigPanSNP_Fst/12WLR_18WLR.fst_0.1.ID",outpre="SNP"  06.SelectionPrefernece_Fst.sbatch

# usage:sbatch --export=input="./LR_Gen",fst="WLR_indel01freq.txt.frq.strat",delta1="0.46",delta2="0.44",id1="./PigPanINDEL_Fst/03WLR_12WLR.fst_0.1.ID",id2="./PigPanIndel_Fst/12WLR_18WLR.fst_0.1.ID",outpre="INDEL"  06.SelectionPrefernece_Fst.sbatch

# Job name:
#SBATCH --job-name=06.SelectionPrefernece_Fst
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-low-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7000MB
#SBATCH --output=logs/06.SelectionPrefernece_Fst_%j.out
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

echo  $(date +"%Y-%m-%d %T") start ${outpre}
cd ${input}

Rscript ${script}/SelectionPrefernece_Fst.R ${fst} ${delta1} ${delta2} ${id1} ${id2} ${outpre} INDEL_large30bp.id

echo  $(date +"%Y-%m-%d %T") all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID