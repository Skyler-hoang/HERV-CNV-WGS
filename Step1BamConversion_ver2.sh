#!/bin/bash

# Script to convert all BAM files from hg19 to hg38 using CrossMap
# Author: Claude
# Date: March 21, 2025

# Define paths
INPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder"
OUTPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/output/hg38_bams"
CHAIN_FILE="/home/hhh38/herv_analysis/reference/hg19ToHg38.over.chain"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Load modules (adjust based on your HPC environment)
# Uncomment and modify these lines as needed for your HPC
# module load crossmap
# module load samtools


THREADS=$(nproc)
MEMORY=64


module load bzip2
module load CrossMap
module load Python/3.9.5-GCCcore-10.3.0

# Log file
LOG_FILE="$OUTPUT_DIR/conversion_log_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Initialize log
log_message "Starting BAM conversion from hg19 to hg38"
log_message "Input directory: $INPUT_DIR"
log_message "Output directory: $OUTPUT_DIR"
log_message "Chain file: $CHAIN_FILE"

# Check if the input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    log_message "ERROR: Input directory does not exist. Exiting."
    exit 1
fi

# Check if the chain file exists
if [ ! -f "$CHAIN_FILE" ]; then
    log_message "ERROR: Chain file does not exist. Exiting."
    exit 1
fi

# Count the number of BAM files
BAM_COUNT=$(find "$INPUT_DIR" -name "*.bam" | wc -l)
log_message "Found $BAM_COUNT BAM files to convert"

# Process each BAM file
COUNTER=0
for BAM_FILE in "$INPUT_DIR"/*.bam; do
    if [ -f "$BAM_FILE" ]; then
        COUNTER=$((COUNTER + 1))
        FILENAME=$(basename "$BAM_FILE")
        OUTPUT_BAM="$OUTPUT_DIR/${FILENAME%.bam}.hg38"
        
        log_message "[$COUNTER/$BAM_COUNT] Converting: $FILENAME"
        log_message "Output will be: $(basename "$OUTPUT_BAM")"
        
        # Run CrossMap
        CrossMap bam -a "$CHAIN_FILE" "$BAM_FILE" "$OUTPUT_BAM"
        
        # Check if conversion was successful
        if [ $? -eq 0 ]; then
            log_message "Conversion successful for $FILENAME"
            
            # Index the output BAM file
            log_message "Indexing $OUTPUT_BAM"
            samtools index "$OUTPUT_BAM"
        else
            log_message "ERROR: Conversion failed for $FILENAME"
        fi
    fi
done

log_message "Conversion complete. $COUNTER files processed."
log_message "Output files are in: $OUTPUT_DIR"

# Print summary
echo "============================================"
echo "Conversion Summary"
echo "============================================"
echo "Total files processed: $COUNTER"
echo "Successful conversions: $(find "$OUTPUT_DIR" -name "*.hg38.bam" | wc -l)"
echo "Check the log file for details: $LOG_FILE"
echo "============================================"
