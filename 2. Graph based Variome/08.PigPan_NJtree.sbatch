#!/bin/bash
# usage:sbatch 08.PigPan_NJtree.sbatch
# Job name:
#SBATCH --job-name=08.PigPan_NJtree
#SBATCH --account=szzx
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-high-mem
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=15000MB
#SBATCH -t 48:00:00
#SBATCH --output=logs/08.PigPan_NJtree_%j.out
### software source
## VCF2Dis: https://github.com/BGI-shenzhen/VCF2Dis
## fneighbor: ftp://emboss.open-bio.org/pub/EMBOSS/ #若edge打不开,可用 filezilla 站点：emboss.open-bio.org,匿名
## consense: http://evolution.genetics.washington.edu/phylip.html
## INSTALL_GUIDE (chinese): https://github.com/hewm2008/VCF2Dis/blob/main/Install.NJ.cn.md
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
conda activate tools_env

input=./
output=./
mkdir -p $output
####SNP
echo  $(date +"%Y-%m-%d %T") filter PanSNP maxmissing 0.9 maf 0.05
bcftools filter -e "F_MISSING > 0.1 ||  MAF < 0.05" ${input}/PanSNP_autosome.vcf.gz --threads 64 -O z -o ${output}/PanSNP_autosome_QC.vcf.gz
tabix -p vcf ${output}/${output}/PanSNP_autosome_QC.vcf.gz

pop=PanSNP_autosome_QC
for (( bootstrap = 0; bootstrap < 100; bootstrap++ )); do

	VCF2Dis -InPut ${output}/${pop}.vcf.gz -OutPut ${pop}.${bootstrap}.mat -Rand 0.25

    fneighbor -datafile ${pop}.${bootstrap}.mat -outfile ${pop}.${bootstrap}.txt -matrixtype s -treetype n -outtreefile ${pop}.${bootstrap}.tre

    echo  $(date +"%Y-%m-%d %T") done ${bootstrap}
    cat ${pop}.*.tre > ALLtree_merge.tre
    fconsense -intreefile ALLtree_merge.tre -outfile out -treeprint Y

    echo  $(date +"%Y-%m-%d %T") construct boostrap nj-tree
    perl  percentageboostrapTree.pl alltree_merge.treefile 100 Final_${pop}_boostrap.tre
done

####Indel
echo  $(date +"%Y-%m-%d %T") filter PanIndel maxmissing 0.9 maf 0.05
bcftools filter -e "F_MISSING > 0.1 ||  MAF < 0.05" ${input}/PanIndel_autosome.vcf.gz --threads 64 -O z -o ${output}/PanIndel_autosome_QC.vcf.gz
tabix -p vcf ${output}/${output}/PanIndel_autosome_QC.vcf.gz

pop=PanIndel_autosome_QC
for (( bootstrap = 0; bootstrap < 100; bootstrap++ )); do

	VCF2Dis -InPut ${output}/${pop}.vcf.gz -OutPut ${pop}.${bootstrap}.mat -Rand 0.25

    fneighbor -datafile ${pop}.${bootstrap}.mat -outfile ${pop}.${bootstrap}.txt -matrixtype s -treetype n -outtreefile ${pop}.${bootstrap}.tre

    echo  $(date +"%Y-%m-%d %T") done ${bootstrap}
    cat ${pop}.*.tre > ALLtree_merge.tre
    fconsense -intreefile ALLtree_merge.tre -outfile out -treeprint Y

    echo  $(date +"%Y-%m-%d %T") construct boostrap nj-tree
    perl  percentageboostrapTree.pl alltree_merge.treefile 100 Final_${pop}_boostrap.tre
done

####SV
echo  $(date +"%Y-%m-%d %T") filter PanSV maxmissing 0.9 maf 0.05
bcftools filter -e "F_MISSING > 0.1 ||  MAF < 0.05" ${input}/PanSV_autosome.vcf.gz --threads 64 -O z -o ${output}/PanSV_autosome_QC.vcf.gz
tabix -p vcf ${output}/${output}/PanSV_autosome_QC.vcf.gz

pop=PanSV_autosome_QC
for (( bootstrap = 0; bootstrap < 100; bootstrap++ )); do

	VCF2Dis -InPut ${output}/${pop}.vcf.gz -OutPut ${pop}.${bootstrap}.mat -Rand 0.25

    fneighbor -datafile ${pop}.${bootstrap}.mat -outfile ${pop}.${bootstrap}.txt -matrixtype s -treetype n -outtreefile ${pop}.${bootstrap}.tre

    echo  $(date +"%Y-%m-%d %T") done ${bootstrap}
    cat ${pop}.*.tre > ALLtree_merge.tre
    fconsense -intreefile ALLtree_merge.tre -outfile out -treeprint Y

    echo  $(date +"%Y-%m-%d %T") construct boostrap nj-tree
    perl  percentageboostrapTree.pl alltree_merge.treefile 100 Final_${pop}_boostrap.tre
done

echo  $(date +"%Y-%m-%d %T")  done

echo "++++++++++++++++++++++++++++++++++++++++"
echo "processs will sleep 30s"
sleep 30
echo "process end at : "
date
rm -rf $CURDIR/nodelist.$SLURM_JOB_ID
