#!/bin/bash

# ERVcaller Analysis Pipeline for Case Western Reserve University HPC
# This script performs ERVcaller analysis for HERV insertions
# It assumes BAM files have already been converted from hg19 to hg38

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# ===================== CONFIGURATION =====================
# Core Directories
WORK_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis"
INPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons"
OUTPUT_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/output"
TEMP_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/temp"
LOG_DIR="/mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/logs"

# Reference Genome Paths
HG38_REFERENCE="$HOME/herv_analysis/reference/hg38/hg38.analysisSet.fa"

# Tool Directories
ERVCALLER_DIR="$HOME/herv_analysis/tools/ERVcaller"
GEVE_FA="$HOME/herv_analysis/reference/gEVE.fa"

# Computational Resources
THREADS=$(nproc)  # Dynamically detect available cores

# ERVcaller Parameters
READ_LENGTH=150
INSERT_SIZE=500
INSERT_STD=100
SPLIT_READ_LEN=20

# ===================== LOGGING SETUP =====================
# Ensure log directory exists
mkdir -p "$LOG_DIR"

# Create timestamped log file
LOG_FILE="${LOG_DIR}/ervcaller_analysis_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to log file and terminal
exec > >(tee -a "$LOG_FILE") 2>&1

# ===================== PREREQUISITE CHECKS =====================
# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check file existence
file_exists() {
    [ -f "$1" ]
}

# Function to check directory existence
dir_exists() {
    [ -d "$1" ]
}

# Validate environment
validate_environment() {
    local errors=0

    # Check required tools
    for tool in perl samtools bwa; do
        if ! command_exists "$tool"; then
            echo "ERROR: $tool is not installed or not in PATH"
            ((errors++))
        fi
    done

    # Check required directories
    for dir in "$WORK_DIR" "$OUTPUT_DIR" "$ERVCALLER_DIR"; do
        if ! dir_exists "$dir"; then
            echo "ERROR: Directory $dir does not exist"
            ((errors++))
        fi
    done

    # Check required files
    for file in "$HG38_REFERENCE" "$GEVE_FA"; do
        if ! file_exists "$file"; then
            echo "ERROR: File $file does not exist"
            ((errors++))
        fi
    done

    # Check if ERVcaller main script exists
    if ! file_exists "${ERVCALLER_DIR}/ERVcaller_v1.4.pl"; then
        echo "ERROR: ERVcaller script not found at ${ERVCALLER_DIR}/ERVcaller_v1.4.pl"
        ((errors++))
    fi

    # Check if hg38 BAM files directory exists
    if ! dir_exists "${OUTPUT_DIR}/hg38_bams"; then
        echo "ERROR: hg38 BAM files directory not found at ${OUTPUT_DIR}/hg38_bams"
        ((errors++))
    fi

    return $errors
}

# ===================== DIRECTORY PREPARATION =====================
prepare_directories() {
    mkdir -p "${OUTPUT_DIR}/ervcaller_analysis"
    mkdir -p "${TEMP_DIR}/ervcaller"
}

# ===================== MAIN FUNCTION =====================
run_ervcaller_analysis() {
    echo "========================================================"
    echo "ERVcaller Analysis - Starting at $(date)"
    echo "========================================================"

    # Validate environment
    if ! validate_environment; then
        echo "ERROR: Environment validation failed. Exiting."
        exit 1
    fi

    # Prepare directories
    prepare_directories

    # Set up ERVcaller environment
    export PATH=$PATH:${ERVCALLER_DIR}
    export PATH=$PATH:${ERVCALLER_DIR}/Scripts

    # Create a file listing all BAM files
    echo "Creating BAM list for ERVcaller..."
    find "${OUTPUT_DIR}/hg38_bams" -name "*.bam" > "${WORK_DIR}/bam_list.txt"

    # Count total files
    TOTAL_FILES=$(wc -l < "${WORK_DIR}/bam_list.txt")
    echo "Found $TOTAL_FILES BAM files to process"

    # Run ERVcaller
    cd "$ERVCALLER_DIR"
    echo "Running ERVcaller..."
    perl ${ERVCALLER_DIR}/ERVcaller_v1.4.pl \
        -i HERV_analysis \
        -f .list \
        -H "$HG38_REFERENCE" \
        -T "$GEVE_FA" \
        -I "${WORK_DIR}" \
        -O "${OUTPUT_DIR}/ervcaller_analysis" \
        -t "$THREADS" \
        -S "$SPLIT_READ_LEN" \
        -BWA_MEM \
        -G \
        -m \
        -i "${WORK_DIR}/bam_list.txt"

    # Check if ERVcaller completed successfully
    if [ $? -ne 0 ]; then
        echo "ERROR: ERVcaller failed to complete successfully"
        exit 1
    fi

    # Combine results if there are multiple samples
    if [ $TOTAL_FILES -gt 1 ]; then
        echo "Combining ERVcaller results from multiple samples..."
        
        # Create a sample list for combining VCFs
        ls "${OUTPUT_DIR}/ervcaller_analysis/"*".vcf" > "${WORK_DIR}/sample_vcf_list.txt"
        
        # Combine VCF files
        perl ${ERVCALLER_DIR}/Scripts/Combine_VCF_files.pl \
            -l "${WORK_DIR}/sample_vcf_list.txt" \
            -o "${OUTPUT_DIR}/ervcaller_analysis/combined.vcf"
        
        # Calculate non-TE insertions for each sample
        find "${OUTPUT_DIR}/hg38_bams" -name "*.bam" | while read bam_file; do
            sample=$(basename "$bam_file" .hg38.bam)
            
            echo "Processing non-TE calculations for sample: $sample"
            perl ${ERVCALLER_DIR}/Scripts/Calculate_reads_among_nonTE_locations.pl \
                -i "${OUTPUT_DIR}/ervcaller_analysis/combined.vcf" \
                -S "$sample" \
                -o "${OUTPUT_DIR}/ervcaller_analysis/${sample}.nonTE" \
                -b "$bam_file" \
                -s paired-end \
                -l "$INSERT_SIZE" \
                -L "$INSERT_STD" \
                -r "$READ_LENGTH" \
                -t "$THREADS"
        done
        
        # Combine non-TE results
        cat ${OUTPUT_DIR}/ervcaller_analysis/*.nonTE > "${OUTPUT_DIR}/ervcaller_analysis/nonTE_allsamples"
        
        # Generate final genotypes
        perl ${ERVCALLER_DIR}/Scripts/Distinguish_nonTE_from_missing_genotype.pl \
            -n "${OUTPUT_DIR}/ervcaller_analysis/nonTE_allsamples" \
            -v "${OUTPUT_DIR}/ervcaller_analysis/combined.vcf" \
            -o "${OUTPUT_DIR}/ervcaller_analysis/final-results.vcf"
    else
        # For single sample, just copy the VCF
        echo "Single sample detected, copying VCF to final-results.vcf"
        cp "${OUTPUT_DIR}/ervcaller_analysis/"*".vcf" "${OUTPUT_DIR}/ervcaller_analysis/final-results.vcf"
    fi

    # Convert VCF to tab-delimited format for easier analysis
    echo "Converting VCF to tab-delimited format..."
    grep -v "^#" "${OUTPUT_DIR}/ervcaller_analysis/final-results.vcf" | \
        awk -F'\t' 'BEGIN {OFS="\t"; print "Chromosome\tPosition\tType\tSample\tLength\tStatus\tGenotype\tGQ"} 
        {
            split($8, info, ";");
            for (i in info) {
                if (info[i] ~ /^INFOR/) {
                    split(info[i], infor_parts, "=");
                    split(infor_parts[2], te_info, ",");
                    type = te_info[1];
                    length = te_info[3];
                    status = te_info[6];
                }
                if (info[i] ~ /^GQ/) {
                    split(info[i], gq_parts, "=");
                    gq = gq_parts[2];
                }
            }
            
            for (i=10; i<=NF; i++) {
                split($i, gt, ":");
                if (gt[1] != "./." && gt[1] != "0/0") {
                    print $1, $2, type, $i, length, status, gt[1], gq;
                }
            }
        }' > "${OUTPUT_DIR}/ervcaller_analysis/final-results.txt"

    # Generate summary statistics
    echo "Generating summary statistics..."
    echo -e "Sample\tTotal_Insertions\tHERV_Insertions\tOther_Insertions" > "${OUTPUT_DIR}/ervcaller_analysis/summary_stats.txt"
    
    awk -F'\t' 'NR>1 {samples[$4]++; if($3 ~ /HERV/) herv[$4]++; else other[$4]++} 
    END {
        for(s in samples) {
            h = herv[s] ? herv[s] : 0;
            o = other[s] ? other[s] : 0;
            print s"\t"samples[s]"\t"h"\t"o
        }
    }' "${OUTPUT_DIR}/ervcaller_analysis/final-results.txt" >> "${OUTPUT_DIR}/ervcaller_analysis/summary_stats.txt"

    echo "========================================================"
    echo "ERVcaller Analysis - Completed at $(date)"
    echo "========================================================"
    echo "Results available at: ${OUTPUT_DIR}/ervcaller_analysis/"
    echo "Summary statistics: ${OUTPUT_DIR}/ervcaller_analysis/summary_stats.txt"
    echo "========================================================"
}

# Run the main function
run_ervcaller_analysis
