#!/bin/bash
#SBATCH --exclusive
# usage:sbatch 7.vg_autoindex_0.8.sbatch
# Job name:
#SBATCH --job-name=7.vg_autoindex_0.8.sbatch
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=hebhcnormal01
###############SBATCH --mem=1000000M
#SBATCH --output=7.vg_autoindex_0.8_%j.out

CONSTRUCTION="./"
LARGE_DISK_TMP=${CONSTRUCTION}/7.vg_autoindex_0.8_tmp
mkdir -p ${LARGE_DISK_TMP}
REFERENCE=W64_addSscrofa11.1ChrYChrMT_removesyncContig.primary_assembly.fa
GRAPHS=${CONSTRUCTION}/7.vg_autoindex_0.8

echo $(date +"%Y-%m-%d %T") phase
zcat ${CONSTRUCTION}/6.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8.vcf.gz | awk -v OFS="\t" '
    BEGIN {FS=OFS; srand()}
    {
        if (/^#/) {
            print;
        } else {
            # 处理基因型信息
            for (i=10; i<=NF; i++) {
                gsub(/\//, "|", $i);
            }
            print; # 打印处理后的行
        }
    }' > ${CONSTRUCTION}/7.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8_transphased.vcf
bgzip -f --threads 24 ${CONSTRUCTION}/7.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8_transphased.vcf
tabix -f -p vcf ${CONSTRUCTION}/7.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8_transphased.vcf.gz
echo $(date +"%Y-%m-%d %T") phase done

mkdir -p ${GRAPHS}
# Build the VG graphs
echo $(date +"%Y-%m-%d %T") Executing vg autoindex
/bin/time -v -p ${CONSTRUCTION}/vg autoindex --workflow giraffe -r ${REFERENCE} -v ${CONSTRUCTION}/7.ONT_Chromosome_Scaffold_Contig_merge_RemoveOverlapSV_all_0.8_transphased.vcf.gz -T ${LARGE_DISK_TMP} -t 64 -V 2 -p ${GRAPHS}/out
echo $(date +"%Y-%m-%d %T") done vg auto

echo $(date +"%Y-%m-%d %T") Executing vg snarl
/bin/time -v -p ${CONSTRUCTION}/vg snarls --threads 64 ${GRAPHS}/out.giraffe.gbz > ${GRAPHS}/out.snarls
echo $(date +"%Y-%m-%d %T") done vg snarl

