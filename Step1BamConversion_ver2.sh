#!/bin/bash

# HERV Analysis Pipeline for Case Western Reserve University HPC
# This script performs:
# 1) BAM conversion from hg19 to hg38
# 2) GenomeSTRIP analysis for duplications & deletions
# 3) ERVcaller analysis for insertions
# 4) ANNOVAR annotation with gEVE
# 5) Integration of all HERV variants

set -e  # Exit on error

# ===================== CONFIGURATION =====================
# Paths - Customized for CWRU HPC
WORK_DIR="$HOME/herv_analysis"
INPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons"  # Original BAM location
OUTPUT_DIR="$WORK_DIR/output"
TEMP_DIR="$WORK_DIR/temp"
LOG_DIR="$WORK_DIR/logs"

# Reference genomes - Update these paths when located
HG19_REFERENCE="/mnt/pan/SOM_EPBI_SKI/Reference/Hg19.fa"
HG38_REFERENCE="$WORK_DIR/reference/hg38.fa"
CHAIN_FILE="$WORK_DIR/reference/hg19ToHg38.over.chain"

# Tool directories - Will need to set up in home directory
GENOMESTRIP_DIR="$WORK_DIR/tools/genomestrip"
ANNOVAR_DIR="$WORK_DIR/tools/annovar"
ERVCALLER_DIR="$WORK_DIR/tools/ERVcaller"
GEVE_GTF="$WORK_DIR/reference/gEVE.gtf"
GEVE_FA="$WORK_DIR/reference/gEVE.fa"

# Computational resources - Adjust based on HPC allocation
THREADS=8
MAX_MEMORY="16g"

# ERVcaller parameters
READ_LENGTH=150
INSERT_SIZE=500
INSERT_STD=100
SPLIT_READ_LEN=20

# GenomeSTRIP parameters
WINDOW_SIZE=1000
WINDOW_OVERLAP=500
REF_GAP_LENGTH=1000
BOUNDARY_PRECISION=100
MIN_REFINED_LENGTH=500

# ===================== CREATE DIRECTORIES =====================
mkdir -p "${OUTPUT_DIR}/hg38_bams"
mkdir -p "${OUTPUT_DIR}/genomestrip/preprocess"
mkdir -p "${OUTPUT_DIR}/genomestrip/sv_discovery"
mkdir -p "${OUTPUT_DIR}/genomestrip/cnv_discovery"
mkdir -p "${OUTPUT_DIR}/genomestrip/logs"
mkdir -p "${OUTPUT_DIR}/annovar"
mkdir -p "${OUTPUT_DIR}/ervcaller_analysis"
mkdir -p "${OUTPUT_DIR}/results"
mkdir -p "${TEMP_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "$WORK_DIR/reference"
mkdir -p "$WORK_DIR/tools"

# Log file
LOG_FILE="${LOG_DIR}/herv_analysis_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================"
echo "HERV Analysis Pipeline - Starting at $(date)"
echo "========================================================"
echo "Working directory: ${WORK_DIR}"
echo "Input directory: ${INPUT_DIR}"
echo "Number of threads: ${THREADS}"
echo "Maximum memory: ${MAX_MEMORY}"

# ===================== LOAD REQUIRED MODULES =====================
echo "Loading required modules..."
module load Python/3.9.5-GCCcore-10.3.0
module load bzip2
echo "Modules loaded successfully."

# ===================== CHECK PREREQUISITES =====================
echo "Checking for necessary tools and reference files..."

# Check for required tools
if ! command -v CrossMap &> /dev/null; then
    echo "CrossMap not found. Installing..."
    pip install --user CrossMap
fi

if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found. Please install samtools before proceeding."
    exit 1
fi

# Check for required reference files
if [ ! -f "$HG38_REFERENCE" ]; then
    echo "ERROR: hg38 reference genome not found at $HG38_REFERENCE"
    echo "Please download it before proceeding."
    exit 1
fi

if [ ! -f "$CHAIN_FILE" ]; then
    echo "ERROR: Chain file not found at $CHAIN_FILE"
    echo "Please download it before proceeding."
    exit 1
fi

echo "All prerequisites checked."

# ===================== STEP 1: BAM CONVERSION =====================
echo "STEP 1: Converting BAM files from hg19 to hg38"

# Get all BAM files
BAM_FILES=($(find "$INPUT_DIR" -maxdepth 1 -name "*.bam"))
TOTAL_FILES=${#BAM_FILES[@]}
echo "Found $TOTAL_FILES BAM files to process"

# Function to convert a single BAM file
convert_bam_to_hg38() {
    local input_bam="$1"
    local filename=$(basename "$input_bam" .bam)
    local output_bam="${OUTPUT_DIR}/hg38_bams/${filename}.hg38.bam"
    local temp_dir="${TEMP_DIR}/${filename}_$$"
    
    echo "Processing $filename"
    
    # Skip if output already exists
    if [ -f "$output_bam" ] && [ -f "${output_bam}.bai" ]; then
        echo "Skipping conversion - output file already exists"
        return 0
    fi
    
    # Create temporary directory
    mkdir -p "$temp_dir"
    
    # Extract sample name from BAM header
    local sample_name=$(samtools view -H "$input_bam" | grep "^@RG" | sed "s/.*SM:\([^ \t]*\).*/\1/" | head -1)
    sample_name=${sample_name:-$filename}
    
    echo "Extracting reads to SAM format..."
    # Process in a pipeline without creating intermediate files where possible
    if samtools view -h "$input_bam" > "${temp_dir}/${filename}.sam"; then
        echo "Converting coordinates to hg38..."
        if CrossMap sam "$CHAIN_FILE" "${temp_dir}/${filename}.sam" "$HG38_REFERENCE" "${temp_dir}/${filename}.hg38.sam"; then
            echo "Converting to BAM format and sorting..."
            if samtools view -bS "${temp_dir}/${filename}.hg38.sam" | \
               samtools sort -@ "$THREADS" -o "${temp_dir}/${filename}.hg38.sorted.bam"; then
                echo "Indexing BAM file..."
                if samtools index "${temp_dir}/${filename}.hg38.sorted.bam"; then
                    echo "Fixing read groups..."
                    if samtools addreplacerg -r "ID:${sample_name}" -r "SM:${sample_name}" \
                       -o "$output_bam" "${temp_dir}/${filename}.hg38.sorted.bam" && \
                       samtools index "$output_bam"; then
                        echo "Conversion completed successfully"
                        # Remove intermediate files to save space
                        rm -f "${temp_dir}/${filename}.sam" "${temp_dir}/${filename}.hg38.sam"
                        rm -f "${temp_dir}/${filename}.hg38.sorted.bam" "${temp_dir}/${filename}.hg38.sorted.bam.bai"
                        rmdir "$temp_dir"
                        return 0
                    else
                        echo "ERROR: Failed to add read groups or index final BAM"
                    fi
                else
                    echo "ERROR: Failed to index sorted BAM"
                fi
            else
                echo "ERROR: Failed to convert to BAM or sort"
            fi
        else
            echo "ERROR: Failed to convert coordinates to hg38"
        fi
    else
        echo "ERROR: Failed to extract reads to SAM format"
    fi
    
    # If we reach here, an error occurred
    echo "Conversion failed for $filename"
    return 1
}

# Process each BAM file one at a time
for ((i=0; i<TOTAL_FILES; i++)); do
    echo "[$((i+1))/$TOTAL_FILES] Processing ${BAM_FILES[$i]}"
    convert_bam_to_hg38 "${BAM_FILES[$i]}"
    echo "------------------------------------------------------"
done

echo "STEP 1 completed: BAM conversion"
