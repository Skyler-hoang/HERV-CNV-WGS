#$ perl user_installed_path/ERVcaller_v.1.4.pl -i TE_seq -f .bam -H hg38.fa -T TE_consensus.fa -r 150 -I folder_of_input_data -O folder_for_output_files -t 12 -S 50 -BWA_MEM -G

#!/bin/bash

# Script to perform ERVcaller
# Author: 
# Date: March 21, 2025

# Define paths
INPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams"

OUTPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/output/ervcaller_analysis"

CONSENSUS_FILE="/$HOME/hhh38/herv_analysis/tools/ERVcaller/ERVcaller-1.4/TE_consensus.fa"

ERV_CALLER="/$HOME/hhh38/herv_analysis/tools/ERVcaller/ERVcaller-1.4/ERVcaller_v1.4.pl"




# Log file
LOG_FILE="$OUTPUT_DIR/conversion_log_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Initialize log
log_message "Starting ERVcaller"
log_message "Input directory: $INPUT_DIR"
log_message "Output directory: $OUTPUT_DIR"
log_message "Consensus file: $CONSENSUS_FILE"
log_message "ERVcaller .pl file": $CONSENSUS_FILE"



# Process each BAM file
COUNTER=0
for BAM_FILE in "$INPUT_DIR"/*.bam; do
    if [ -f "$BAM_FILE" ]; then
        COUNTER=$((COUNTER + 1))
        FILENAME=$(basename "$BAM_FILE")
        OUTPUT_BAM="$OUTPUT_DIR/${FILENAME%.bam}.hg38.bam"
        
        log_message "[$COUNTER/$BAM_COUNT] Converting: $FILENAME"
        log_message "Output will be: $(basename "$OUTPUT_BAM")"
        
        # Run ERVcaller
        perl $ERV_CALLER -i $BAM_FILE -f .bam -H hg38.fa -T $CONSENSUS_FILE -r 150 -I $INPUT_DIR -O $OUTPUT_DIR -t 16 -S 50 -BWA_MEM -G        
    fi
done

log_message "ERVCaller analysis complete. $COUNTER files processed."
log_message "Output files are in: $OUTPUT_DIR"

# Print summary
echo "============================================"
echo "Conversion Summary"
echo "============================================"
echo "Total files processed: $COUNTER"
echo "Successful analyses: $(find "$OUTPUT_DIR" -name "*.hg38.bam" | wc -l)"
echo "Check the log file for details: $LOG_FILE"
echo "============================================"
