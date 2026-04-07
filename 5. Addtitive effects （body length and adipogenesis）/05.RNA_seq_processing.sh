# ==============================================================================
# RNA-Seq Processing Pipeline
# Performs quality control, alignment, and quantification using fastp, STAR, and featureCounts
# ==============================================================================

# --------------------------
# Section 1: Configuration
# --------------------------

# Input/Output paths
METADATA_FILE="metadata.txt"          # Sample metadata file
INFQ_PATH="raw_data"                  # Input FASTQ directory
OUT_PATH="processed_data"             # Output directory
STAR_INDEX="1.ref_GTF_index"          # STAR genome index directory
REF_FASTA="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"  # Reference genome
GTF_FILE="Sus_scrofa.Sscrofa11.1.110.gtf"           # Gene annotation

# --------------------------
# Section 2: Sample Processing Setup
# --------------------------

# Read sample IDs from metadata file
SAMPLES=($(awk '{print $1}' ${METADATA_FILE}))
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample: ${SAMPLE_ID}"

# Define file paths
INPUT_FQ1="${INFQ_PATH}/${SAMPLE_ID}_1.fq.gz"
INPUT_FQ2="${INFQ_PATH}/${SAMPLE_ID}_2.fq.gz"
CLEAN_FQ1="${OUT_PATH}/${SAMPLE_ID}_1.clean.fq.gz"
CLEAN_FQ2="${OUT_PATH}/${SAMPLE_ID}_2.clean.fq.gz"
HTML_REPORT="${OUT_PATH}/output_html/${SAMPLE_ID}.html"
JSON_REPORT="${OUT_PATH}/output_html/${SAMPLE_ID}.json"

# --------------------------
# Section 3: Quality Control with fastp
# --------------------------

if [ ! -s ${JSON_REPORT} ] || [ ! -s ${CLEAN_FQ1} ] || [ ! -s ${CLEAN_FQ2} ]; then
    echo "Running fastp quality control..."
    
    fastp \
        -i ${INPUT_FQ1} -I ${INPUT_FQ2} \
        -o ${CLEAN_FQ1} -O ${CLEAN_FQ2} \
        -h ${HTML_REPORT} -j ${JSON_REPORT} \
        --detect_adapter_for_pe \
        -w ${SLURM_CPUS_PER_TASK} \
        --fix_mgi_id \
        --correction

    # Check fastp exit status
    if [ $? -eq 0 ]; then
        echo "++++++++++++++++++++++++++++++++++++++++++"
        echo "$(date) Successfully processed ${SAMPLE_ID} with fastp"
        echo "++++++++++++++++++++++++++++++++++++++++++"
    else
        echo "++++++++++++++++++++++++++++++++++++++++++"
        echo "$(date) ERROR processing ${SAMPLE_ID} with fastp"
        echo "++++++++++++++++++++++++++++++++++++++++++"
        exit 1
    fi
else
    echo "++++++++++++++++++++++++++++++++++++++++++"
    echo "$(date) ${SAMPLE_ID} already processed - skipping"
    echo "++++++++++++++++++++++++++++++++++++++++++"
fi

# --------------------------
# Section 4: Genome Alignment with STAR
# --------------------------

# Generate genome index if not exists
if [ ! -d "${STAR_INDEX}" ]; then
    echo "Creating STAR genome index..."
    time STAR \
        --runThreadN ${SLURM_CPUS_PER_TASK} \
        --runMode genomeGenerate \
        --genomeDir ${STAR_INDEX} \
        --genomeFastaFiles ${REF_FASTA} \
        --sjdbGTFfile ${GTF_FILE} \
        --sjdbOverhang 149
fi

# Align reads
echo "Aligning ${SAMPLE_ID} with STAR..."
STAR \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --genomeDir ${STAR_INDEX} \
    --sjdbGTFfile ${GTF_FILE} \
    --twopassMode Basic \
    --quantMode TranscriptomeSAM \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 20 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMismatchNmax 999 \
    --outFilterMatchNmin 0 \
    --outTmpDir tmp/${SAMPLE_ID} \
    --readFilesIn ${CLEAN_FQ1} ${CLEAN_FQ2} \
    --outFileNamePrefix ${SAMPLE_ID}_

# --------------------------
# Section 5: featureCounts
# --------------------------

echo "Quantifying gene expression for ${SAMPLE_ID}..."
time featureCounts \
    -p \
    -C \
    -s 2 \
    -t exon \
    -g gene_id \
    -T ${SLURM_CPUS_PER_TASK} \
    -a ${GTF_FILE} \
    -o ${OUT_PATH}/${SAMPLE_ID}_counts.txt \
    ${SAMPLE_ID}.sorted.bam

echo "Process finished for ${SAMPLE_ID}"