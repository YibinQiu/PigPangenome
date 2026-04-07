#!/bin/bash
# usage: sbatch hiblup_grm.sbatch

path1=./50K/1.geno
path2=./50K/2.grm
HIBLUP=./software/hiblup/hiblup

genotype=${path1}/geno
xrm_output=${path2}/50K

# Generate GA,GD files
$HIBLUP --make-xrm --code-method 4 --bfile ${genotype} \
        --add --dom --thread 4 --out $xrm_output

# GRM conversion format
$HIBLUP --trans-xrm --xrm ${xrm_output}.GA --out ${path2}/50K_G
$HIBLUP --trans-xrm --xrm ${xrm_output}.GD --out ${path2}/50K_D

