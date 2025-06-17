#!/bin/bash
# usage: sbatch hiblup_grm.sbatch

path1=./50K/1.geno
path2=./50K/2.grm
HIBLUP=./software/hiblup/hiblup

for i in 1191 996 673 520 517
do
genotype=${path1}/${i}_id
xrm_output=${path2}/${i}_50K

# Generate GA,GD files
$HIBLUP --make-xrm --code-method 4 --bfile ${genotype} \
        --add --dom --thread 4 --out $xrm_output

# GRM conversion format
$HIBLUP --trans-xrm --xrm ${xrm_output}.GA --out ${path2}/${i}_50K_G
$HIBLUP --trans-xrm --xrm ${xrm_output}.GD --out ${path2}/${i}_50K_D
done
