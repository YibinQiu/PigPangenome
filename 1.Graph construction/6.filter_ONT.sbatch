#!/bin/bash
# usage:sbatch 6.filter_ONT.sbatch
# Job name:
#SBATCH --job-name=6.filter_ONT.sbatch
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
#SBATCH --output=6.filter_ONT_%j.out
ref=reference
input=raw_vcf
output=1.ONT_highConfidence

module load anaconda3
module load bcftools
export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C
source activate jasmine

mkdir -p ${output}

while read sample
do
echo $(date +"%Y-%m-%d %T") start ${sample} sensetive filter

bcftools view -i 'FILTER="PASS"' ${input}/${sample}.change_name_add_refseq.vcf | \
bcftools filter -i '((SVTYPE="DEL" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000) || (SVTYPE="DUP" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000) || (SVTYPE="INV" && ABS(SVLEN)>=30 && ABS(SVLEN)<=100000) || (SVTYPE="INS" && ABS(SVLEN)>=30 && ABS(SVLEN)<=20000)) && INFO/SUPPORT>=5' | \
grep -v 'IMPRECISE;' | bcftools sort | awk -v OFS="\t" '
    BEGIN {FS=OFS; srand()}
    {
        if (/^#/) {
            print;
        } else {
            svtype = ""; split($8, info, ";");
            for (i in info) {
                if (info[i] ~ /^SVTYPE=/) {
                    split(info[i], svtype_parts, "=");
                    svtype = svtype_parts[2];
                }
            }
            if (svtype == "INV") {
                $4 = substr($4, 1, 1);
                $5 = "<INV>";
            }
            if (svtype == "DUP") {
                $4 = substr($4, 1, 1);
                $5 = "<DUP>";
            }
            print;
        }
    }' > ${output}/${sample}.highConfidence.vcf

#echo $(date +"%Y-%m-%d %T") start ${sample} dup to ins
#jasmine --comma_filelist file_list=${output}/${sample}.highConfidence.vcf --leave_breakpoints --dup_to_ins genome_file=${ref}/W64_addSscrofa11.1ChrY.primary_assembly.fa --preprocess_only out_dir=${output} threads=24 samtools_path=samtools max_dup_length=100000 out_file=${output}/${sample}.highConfidence_dupToIns.vcf &> ${output}/${sample}.highConfidence_dupToIns.vcf.log

echo $(date +"%Y-%m-%d %T") start ${sample} pre norm
jasmine --preprocess_only --pre_normalize --leave_breakpoints --comma_filelist file_list=${output}/${sample}.highConfidence.vcf out_dir=${output} out_file=${output}/${sample}.highConfidence_normalizeTypes.vcf &> ${output}/${sample}.highConfidence_normalizeTypes.vcf.log

echo $(date +"%Y-%m-%d %T") start ${sample} intra merge
jasmine --comma_filelist file_list=${output}/${sample}.highConfidence_normalizeTypes.vcf --normalize_type --leave_breakpoints threads=24 kd_tree_norm=2 max_dist=100 --allow_intrasample max_dist_linear=0 min_dist=-1 min_seq_id=0 k_jaccard=9 out_dir=${output} out_file=${output}/${sample}.highConfidence_normalizeTypes.ism.vcf > ${output}/${sample}.highConfidence_normalizeTypes.ism.vcf.log

echo $(date +"%Y-%m-%d %T") start ${sample} sort and rm duplicate sv
cat ${output}/${sample}.highConfidence_normalizeTypes.ism.vcf | awk '
    BEGIN { FS="\t"; OFS="\t" }
    {
        if ($1 ~ /^#/) {
           
            print;
        } else {
          
            chrom = $1; pos = $2; info = $8;
            split(info, infos, ";");
            for (i in infos) {
                if (infos[i] ~ /^SVTYPE=/) {
                    split(infos[i], arr, "=");
                    svtype = arr[2];
                }
                if (infos[i] ~ /^SUPPORT=/) {
                    split(infos[i], arr, "=");
                    support = arr[2];
                }
            }
            key = chrom SUBSEP pos SUBSEP svtype;

        
            if (key in sv) {
                if (support > sv[key]["support"]) {
                  
                    sv[key]["line"] = $0;
                    sv[key]["support"] = support;
                }
            } else {
             
                sv[key]["line"] = $0;
                sv[key]["support"] = support;
            }
        }
    }
    END {
     
        for (key in sv) {
            print sv[key]["line"];
        }
    }' | bcftools sort | awk 'BEGIN {FS="\t"; OFS="\t"}
{
    if ($1 ~ /^#/) {
        
        print $0;
    } else {
         
         ref=$4; alt=$5;
         match($8, /SVTYPE=([^;]+)/, arr_svtype);
         svtype = arr_svtype[1];

         
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
}' | bgzip -f > ${output}/${sample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz

tabix -f -p vcf ${output}/${sample}.highConfidence_normalizeTypes.ism.sorted.vcf.gz

echo $(date +"%Y-%m-%d %T") ${sample} sensetive filter done

done < 0.ONT_sample
