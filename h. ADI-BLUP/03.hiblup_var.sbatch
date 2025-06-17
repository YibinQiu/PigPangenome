#!/bin/bash
# usage: sbatch hiblup_var.sbatch

path1=./phe_GS
path2=./50K/2.grm
path3=./50K/3.var
HIBLUP=./software/hiblup/hiblup
gs_array_task=./gs_array_task.txt

# Read trait array
traits=($(cat ${gs_array_task} | awk '{print $1}'))
geno_num=($(cat ${gs_array_task} | awk '{print $2}'))
trait=${traits[$SLURM_ARRAY_TASK_ID]}
i=${geno_num[$SLURM_ARRAY_TASK_ID]}
echo "正在处理性状: ${trait}"
echo "性状的个体数量: ${i}"

# Define specific files
pheno_file=${path1}/${trait}.txt
xrm_output=${path2}/${i}_50K
result_output=${path3}/${trait}_add_dom_epi

# Estimate the variance components
$HIBLUP --single-trait --pheno $pheno_file --pheno-pos 2 \
        --xrm ${xrm_output}.GA,${xrm_output}.GD,${xrm_output}.GA:${xrm_output}.GA \
        --vc-method AI --ai-maxit 1000 --thread 4 --out $result_output

echo "性状 ${trait} 的分析完成，结果保存到 ${result_output}"
