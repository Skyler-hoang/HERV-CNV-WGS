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
# Will need to ask Ricky for these locations
HG19_REFERENCE="/mnt/pan/SOM_EPBI_SKI/Reference/Hg19.fa"
HG38_REFERENCE="$WORK_DIR/reference/hg38.fa"   # Will need to download
CHAIN_FILE="$WORK_DIR/reference/hg19ToHg38.over.chain"  # Will need to download

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

# Check for necessary tools and download reference files
echo "Checking for necessary tools and reference files..."

# Set up Python environment
echo "Setting up Python environment..."
module load Python/3.9.5-GCCcore-10.3.0

# Install required Python packages if needed
pip list | grep -q CrossMap || pip install --user CrossMap
echo "Python environment set up."

# Download hg38 reference and chain file if needed
if [ ! -f "$HG38_REFERENCE" ]; then
    echo "hg38 reference genome not found. Downloading..."
    # This is a placeholder - need to implement actual download
    echo "WARNING: You will need to manually download hg38 reference to $HG38_REFERENCE"
fi

if [ ! -f "$CHAIN_FILE" ]; then
    echo "Chain file not found. Downloading..."
    # This is a placeholder - need to implement actual download
    echo "WARNING: You will need to manually download the chain file to $CHAIN_FILE"
fi

# Check for samtools - try to install locally if not available
if ! command -v samtools &> /dev/null; then
    echo "samtools not found in system path. Will need to install locally."
    # This is a placeholder - need to implement actual install
    echo "WARNING: You will need to manually install samtools locally"
fi

# ===================== STEP 1: BAM CONVERSION =====================
echo "STEP 1: Converting BAM files from hg19 to hg38"
echo "WARNING: This step cannot proceed until samtools is available and reference files are set up."
echo "Please install the required tools and reference files before proceeding."

# This is a placeholder for the conversion step
echo "To continue with the actual conversion, uncomment the code below after installing prerequisites."

: '
# Get all BAM files
BAM_FILES=($(find "$INPUT_DIR" -maxdepth 1 -name "*.bam"))
TOTAL_FILES=${#BAM_FILES[@]}
echo "Found $TOTAL_FILES BAM files to process"

for ((i=0; i<TOTAL_FILES; i++)); do
    input_bam="${BAM_FILES[$i]}"
    filename=$(basename "$input_bam" .bam)
    output_bam="${OUTPUT_DIR}/hg38_bams/${filename}.hg38.bam"
    
    echo "[$((i+1))/$TOTAL_FILES] Processing $filename"
    
    # Skip if output already exists
    if [ -f "$output_bam" ] && [ -f "${output_bam}.bai" ]; then
        echo "Skipping conversion - output file already exists"
        continue
    fi
    
    # Extract sample name from BAM header
    sample_name=$(samtools view -H "$input_bam" | grep "^@RG" | sed "s/.*SM:\([^ \t]*\).*/\1/" | head -1)
    if [ -z "$sample_name" ]; then
        sample_name="$filename"
    fi
    
    # Create temp directory for this file
    file_temp_dir="${TEMP_DIR}/${filename}_$$"
    mkdir -p "$file_temp_dir"
    
    # Extract reads to SAM format
    echo "[$((i+1))/$TOTAL_FILES] Extracting reads to SAM..."
    samtools view -h "$input_bam" > "${file_temp_dir}/${filename}.sam"
    
    # Convert coordinates from hg19 to hg38
    echo "[$((i+1))/$TOTAL_FILES] Converting coordinates to hg38..."
    CrossMap.py sam "$CHAIN_FILE" "${file_temp_dir}/${filename}.sam" "$HG38_REFERENCE" "${file_temp_dir}/${filename}.hg38.sam"
    
    # Convert SAM back to BAM
    echo "[$((i+1))/$TOTAL_FILES] Converting to BAM format..."
    samtools view -bS "${file_temp_dir}/${filename}.hg38.sam" > "${file_temp_dir}/${filename}.hg38.bam"
    
    # Sort BAM file
    echo "[$((i+1))/$TOTAL_FILES] Sorting BAM file..."
    samtools sort -@ "$THREADS" "${file_temp_dir}/${filename}.hg38.bam" -o "${file_temp_dir}/${filename}.hg38.sorted.bam"
    
    # Index BAM file
    echo "[$((i+1))/$TOTAL_FILES] Indexing BAM file..."
    samtools index "${file_temp_dir}/${filename}.hg38.sorted.bam"
    
    # Fix read groups
    echo "[$((i+1))/$TOTAL_FILES] Fixing read groups..."
    samtools addreplacerg -r "ID:${sample_name}" -r "SM:${sample_name}" -o "$output_bam" "${file_temp_dir}/${filename}.hg38.sorted.bam"
    
    # Index final BAM file
    samtools index "$output_bam"
    
    # Clean up
    rm -rf "$file_temp_dir"
    
    echo "[$((i+1))/$TOTAL_FILES] Conversion completed"
done
'

echo "STEP 1 is currently commented out. Please uncomment after setting up prerequisites."
echo "Script has created the directory structure in your home directory at: $WORK_DIR"
echo "Please contact Ricky Chan for assistance with reference files and tools."
