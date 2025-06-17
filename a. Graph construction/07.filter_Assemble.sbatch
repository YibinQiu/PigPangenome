#!/bin/bash
# usage:yhbatch 7.filter_Assemble.sbatch
# Job name:
#SBATCH --job-name=7.filter_Assemble.sbatch
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=rhenv
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=24
# Memory per node:
#SBATCH --output=7.filter_Assemble_%j.out
module load anaconda3
module load bcftools
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C
source activate jasmine
output=1.Assemble_highConfidence
mkdir -p ${output}
path=raw_vcf
while read sampleID
do
    echo $(date +"%Y-%m-%d %T") start ${sampleID} sensetive filter
    zcat ${path}/${sampleID}_svmu.vcf.gz > ${output}/${sampleID}_svmu_nohead.vcf
    cat ${path}/fake_head ${output}/${sampleID}_svmu_nohead.vcf > ${output}/${sampleID}_svmu.vcf
    rm ${output}/${sampleID}_svmu_nohead.vcf 
    specialIDs="Bama Bamei meishan Ninghe NX Sscrofa11.1 USMARC NSME Donna W51hap1 W51hap2 W51 S21hap1 S21hap2 S21 W64hap1 W64hap2"

    
    if [[ " $specialIDs " =~ " $sampleID " ]]; then
       
        bcftoolsFilterCmd="-i '((SVTYPE=\"DEL\" && ABS(SVLEN)>=30 && ABS(SVLEN)<=5000000) || (SVTYPE=\"INS\" && ABS(SVLEN)>=30 && ABS(SVLEN)<=5000000))'"
    else
        
        bcftoolsFilterCmd="-i '((SVTYPE=\"DEL\" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000) || (SVTYPE=\"INS\" && ABS(SVLEN)>=30 && ABS(SVLEN)<=20000))'"
    fi

    
    eval bcftools filter $bcftoolsFilterCmd ${output}/${sampleID}_svmu.vcf | \
    awk -v sampleID="$sampleID" 'BEGIN { FS=OFS="\t" }
    {
        if ($1 ~ /^#/) {
            
            print $0;
        } else {
             
             ref=$4; alt=$5;
             match($8, /SVTYPE=([^;]+)/, arr_svtype);
             svtype = arr_svtype[1];
             $3 = sampleID":"svtype":"$1"_"$2;
             
             ref_len = length(ref);
             alt_len = length(alt);

             
             non_atcg_ref = gsub(/[^ATCGatcg]/, "", ref);
             non_atcg_alt = gsub(/[^ATCGatcg]/, "", alt);

             
             ref_ratio = non_atcg_ref / ref_len;
             alt_ratio = non_atcg_alt / alt_len;
             if (svtype == "INS" || svtype == "DEL") {
             
                 if (ref_ratio > 0.01 || alt_ratio > 0.01) {
                 next;
                 }
             }
             print $0; 

        }
    }'| bcftools sort > ${output}/${sampleID}.highConfidence.vcf
    
    
echo $(date +"%Y-%m-%d %T") start ${sampleID} pre norm
jasmine --preprocess_only --pre_normalize --leave_breakpoints --comma_filelist file_list=${output}/${sampleID}.highConfidence.vcf out_dir=${output} out_file=${output}/${sampleID}.highConfidence_normalizeTypes.vcf &> ${output}/${sampleID}.highConfidence_normalizeTypes.vcf.log

echo $(date +"%Y-%m-%d %T") start ${sampleID} intra merge
jasmine --comma_filelist file_list=${output}/${sampleID}.highConfidence_normalizeTypes.vcf --normalize_type --leave_breakpoints threads=24 kd_tree_norm=2 max_dist=100 --allow_intrasample max_dist_linear=0 min_dist=-1 min_seq_id=0 k_jaccard=9 out_dir=${output} out_file=${output}/${sampleID}.highConfidence_normalizeTypes.ism.vcf > ${output}/${sampleID}.highConfidence_normalizeTypes.ism.vcf.log

cat ${output}/${sampleID}.highConfidence_normalizeTypes.ism.vcf | \
awk 'BEGIN { FS=OFS="\t" }
{
    if ($1 ~ /^#/ || $3 !~ /duplicate/) {
        $4 = toupper($4);
        $5 = toupper($5);
        print $0;
    }
}' | bcftools sort |bgzip -f > ${output}/${sampleID}.highConfidence_normalizeTypes.ism.sorted.vcf.gz

tabix -f -p vcf ${output}/${sampleID}.highConfidence_normalizeTypes.ism.sorted.vcf.gz

done < 0.Assemble_sample

