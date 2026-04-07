#!/bin/bash
#usage: sbatch 02.hypertuning.sbatch
#SBATCH --job-name=02.hypertuning
#SBATCH --qos=nolimit
#SBATCH --account=szzx
#SBATCH --nodes=1
#SBATCH --output=logs/02.hypertuning_%A_%a.out
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=7500MB
#SBATCH --cpus-per-task=5
#SBATCH --partition=cpu-high-mem,cpu-low-mem,cpu-6226
#SBATCH --time=148:00:00

# ============================================================
CURDIR=`pwd`
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID
NODES=`scontrol show hostnames $SLURM_JOB_NODELIST`
for i in $NODES
do
echo "$i:$SLURM_NTASKS_PER_NODE" >> $CURDIR/nodelist.$SLURM_JOB_ID
done
# Load required modules
. ~/miniconda3/bin/activate
conda activate ac

# Define variables
train_id=${SLURM_ARRAY_TASK_ID}
echo "Start train_${train_id}"

base_dir="./"
train_dir="${base_dir}/trains_all/train_${train_id}"
mkdir -p "$train_dir"

# Copy necessary files
cp ${base_dir}/raw_phe.csv "$train_dir"
cp ${base_dir}/file_test.csv "$train_dir"
cp ${base_dir}/modified_* "$train_dir"

# Define copymask levels
copymask_levels=(10 20 30 40 50 60 70 80)

# Train models with different copymask amounts
for copymask_level in "${copymask_levels[@]}"; do
    model_path="${train_dir}/all_train_copymask${copymask_level}.pth"
    python3.7 ~/software/AutoComplete/fit.py \
        "${train_dir}/raw_phe.csv" \
        --id_name ID \
        --copymask_amount 0.${copymask_level} \
        --batch_size 2048 \
        --epochs 500 \
        --lr 0.1 \
        --quality \
        --save_model_path "$model_path"
done

# Evaluate models with imputation
# fit.py and bootstrap_r2_statistic.py can be obtained from https://github.com/sriramlab/AutoComplete
for copymask_level in "${copymask_levels[@]}"; do
    echo "Start copymask${copymask_level}"
    copymask_dir="${train_dir}/copymask${copymask_level}"
    mkdir -p "$copymask_dir"
    cp "${train_dir}/raw_phe.csv" "$copymask_dir"
    cp "${train_dir}/file_test.csv" "$copymask_dir"
    cp ${train_dir}/modified_* "$copymask_dir"
    cp "${train_dir}/all_train_copymask${copymask_level}.pth" "$copymask_dir"
    cp "${train_dir}/all_train_copymask${copymask_level}.pth.json" "$copymask_dir"
    cp "${train_dir}/all_train_copymask${copymask_level}_quality.csv" "$copymask_dir"

    for modified_file in modified_5.csv modified_10.csv modified_20.csv modified_30.csv modified_40.csv modified_50.csv modified_60.csv; do
        echo "####################################"
        imputed_file="${copymask_dir}/${modified_file%.csv}_imputed.csv"
        echo "${copymask_dir}/all_train_copymask${copymask_level}.pth"
        echo "${copymask_dir}/${modified_file}"
        echo "$imputed_file"
        python3.7 ~/software/AutoComplete/fit.py \
            "${copymask_dir}/raw_phe.csv" \
            --id_name ID \
            --impute_using_saved "${copymask_dir}/all_train_copymask${copymask_level}.pth" \
            --impute_data_file "${copymask_dir}/${modified_file}" \
            --output "$imputed_file"
    
        result_file="${copymask_dir}/result_copymask${copymask_level}_${modified_file%.csv}.csv"
        echo "${copymask_dir}/${modified_file}"
        echo "${copymask_dir}/${modified_file%.csv}_imputed.csv"
        echo "$result_file"
        python3.7 ~/software/AutoComplete/bootstrap_r2_statistic.py \
            "${copymask_dir}/file_test.csv" \
            --simulated_data_file "${copymask_dir}/${modified_file}" \
            --imputed_data_file "${copymask_dir}/${modified_file%.csv}_imputed.csv" \
            --num_bootstraps 100 \
            --saveas "$result_file"
    done

done

echo "Completed train_${train_id}"
echo $(date +"%m-%d-%Y-%T"): all done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 60
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID