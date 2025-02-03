#!/bin/bash
set -e
# usage:sbatch -a 1-2598 1.giraffe_genotyping.sbatch 1 2598 totalsample
# Job name:
#SBATCH --job-name=1.giraffe_genotyping.sbatch
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=rhenv
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
# Processors per task:
#SBATCH --cpus-per-task=24
#SBATCH --output=1.giraffe_genotyping_%A_%a_%j.out

start=$1
end=$2
samplefile=$3
sample=(`cat ${samplefile} | awk '{if(NR>='"$start"'&&NR<='"$end"')print}'`)
echo ${sample[SLURM_ARRAY_TASK_ID]}_${SLURM_JOB_ID}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_giraffe
module load samtools

export LANGUAGE="en_US.UTF-8"
export LANG=en_US:zh_CN.UTF-8
export LC_ALL=C
#######################################################################################
LARGE_DISK_TMP="tmp"
export TMPDIR=${LARGE_DISK_TMP}
mkdir -p ${LARGE_DISK_TMP}
input="./"
vgpath="./"

GRAPHS="./"
output="./"
gamfile_store="./"

mkdir -p ${output}
mkdir -p ${gamfile_store}
mkdir -p ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}
echo $(date +"%Y-%m-%d %T") begin copy
cp ${GRAPHS}/out.dist ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.dist
cp ${GRAPHS}/out.giraffe.gbz ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.giraffe.gbz
cp ${GRAPHS}/out.min ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.min
cp ${GRAPHS}/out.snarls ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.snarls
echo $(date +"%Y-%m-%d %T") done copy
################################################################################
echo  $(date +"%Y-%m-%d %T") begin giraffe ${sample[SLURM_ARRAY_TASK_ID]}
/bin/time -v -p ${vgpath}/vg giraffe -p -t 24 -Z ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.giraffe.gbz -d ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.dist -m ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.min \
    -f ${input}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}_1.fq.gz -f ${input}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}_2.fq.gz > ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.gam
echo  $(date +"%Y-%m-%d %T") done giraffe ${sample[SLURM_ARRAY_TASK_ID]}

echo  $(date +"%Y-%m-%d %T") start pack ${sample[SLURM_ARRAY_TASK_ID]}
/bin/time -v -p ${vgpath}/vg pack -t 24 -x ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.giraffe.gbz -g ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.gam -o ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.pack -Q 5
echo  $(date +"%Y-%m-%d %T") done pack ${sample[SLURM_ARRAY_TASK_ID]}

echo  $(date +"%Y-%m-%d %T") start call ${sample[SLURM_ARRAY_TASK_ID]}
/bin/time -v -p ${vgpath}/vg call -t 24 ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.giraffe.gbz -r ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.snarls -k ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.pack -s ${sample[SLURM_ARRAY_TASK_ID]} -a -z > ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf
echo  $(date +"%Y-%m-%d %T") done call ${sample[SLURM_ARRAY_TASK_ID]}

echo $(date +"%Y-%m-%d %T") begin bgzip index ${sample[SLURM_ARRAY_TASK_ID]}.vcf
bgzip -f --threads 24 ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf
tabix -f -p vcf ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf.gz
echo $(date +"%Y-%m-%d %T") done bgzip index ${sample[SLURM_ARRAY_TASK_ID]}.vcf

echo $(date +"%Y-%m-%d %T") executing stats
/bin/time -v -p ${vgpath}/vg stats -a ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.gam >${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.stats
echo  $(date +"%Y-%m-%d %T") done stats ${sample[SLURM_ARRAY_TASK_ID]}

mv ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.stats ${output}/
mv ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf.gz ${output}/
mv ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.vcf.gz.tbi ${output}/
mv ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.pack ${output}/
#######################################################################################
rm ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.dist
rm ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.giraffe.gbz
rm ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.min
rm ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/out.snarls
echo $(date +"%m-%d-%Y-%T"): done remove graph
mv ${TMPDIR}/${sample[SLURM_ARRAY_TASK_ID]}/${sample[SLURM_ARRAY_TASK_ID]}.gam ${gamfile_store}/
echo $(date +"%m-%d-%Y-%T"): all done
