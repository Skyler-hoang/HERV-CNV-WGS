#!/bin/bash

# HERV Analysis Pipeline
# This script performs:
# 1) BAM conversion from hg19 to hg38
# 2) GenomeSTRIP analysis for duplications & deletions
# 3) ERVcaller analysis for insertions
# 4) ANNOVAR annotation with gEVE
# 5) Integration of all HERV variants

set -e  # Exit on error

# ===================== CONFIGURATION =====================
# Paths - MODIFY THESE
WORK_DIR="/path/to/working/directory"
INPUT_DIR="${WORK_DIR}/input"  # Directory with hg19 BAM files
OUTPUT_DIR="${WORK_DIR}/output"
HG19_REFERENCE="/path/to/hg19.fa"
HG38_REFERENCE="/path/to/hg38.fa"
CHAIN_FILE="/path/to/hg19ToHg38.over.chain"
GENOMESTRIP_DIR="/path/to/genomestrip"
ANNOVAR_DIR="/path/to/annovar"
ERVCALLER_DIR="/path/to/ERVcaller_v1.4"
GEVE_GTF="/path/to/gEVE.gtf"
GEVE_FA="/path/to/gEVE.fa"
THREADS=16
MAX_MEMORY="32g"

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

# Log file
LOG_FILE="${WORK_DIR}/herv_analysis.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================"
echo "HERV Analysis Pipeline - Starting at $(date)"
echo "========================================================"

# ===================== STEP 1: BAM CONVERSION =====================
echo "STEP 1: Converting BAM files from hg19 to hg38"

# Get all BAM files
BAM_FILES=($(find "$INPUT_DIR" -name "*.bam"))
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
    sample_name=$(samtools view -H "$input_bam" | grep '^@RG' | sed 's/.*SM:\([^ \t]*\).*/\1/' | head -1)
    if [ -z "$sample_name" ]; then
        sample_name="$filename"
    fi
    
    # Create temp directory for this file
    temp_dir="${WORK_DIR}/temp_${filename}_$$"
    mkdir -p "$temp_dir"
    
    # Extract reads to SAM format
    echo "[$((i+1))/$TOTAL_FILES] Extracting reads to SAM..."
    samtools view -h "$input_bam" > "${temp_dir}/${filename}.sam"
    
    # Convert coordinates from hg19 to hg38
    echo "[$((i+1))/$TOTAL_FILES] Converting coordinates to hg38..."
    CrossMap.py sam "$CHAIN_FILE" "${temp_dir}/${filename}.sam" "$HG38_REFERENCE" "${temp_dir}/${filename}.hg38.sam"
    
    # Convert SAM back to BAM
    echo "[$((i+1))/$TOTAL_FILES] Converting to BAM format..."
    samtools view -bS "${temp_dir}/${filename}.hg38.sam" > "${temp_dir}/${filename}.hg38.bam"
    
    # Sort BAM file
    echo "[$((i+1))/$TOTAL_FILES] Sorting BAM file..."
    samtools sort -@ "$THREADS" "${temp_dir}/${filename}.hg38.bam" -o "${temp_dir}/${filename}.hg38.sorted.bam"
    
    # Index BAM file
    echo "[$((i+1))/$TOTAL_FILES] Indexing BAM file..."
    samtools index "${temp_dir}/${filename}.hg38.sorted.bam"
    
    # Fix read groups
    echo "[$((i+1))/$TOTAL_FILES] Fixing read groups..."
    samtools addreplacerg -r "ID:${sample_name}" -r "SM:${sample_name}" -o "$output_bam" "${temp_dir}/${filename}.hg38.sorted.bam"
    
    # Index final BAM file
    samtools index "$output_bam"
    
    # Clean up
    rm -rf "$temp_dir"
    
    echo "[$((i+1))/$TOTAL_FILES] Conversion completed"
done

echo "STEP 1 completed: All BAM files converted to hg38"

# ===================== STEP 2: GENOMESTRIP ANALYSIS =====================
echo "STEP 2: Running GenomeSTRIP analysis"

# Create necessary files for GenomeSTRIP
echo "Creating metadata and configuration files for GenomeSTRIP..."

# Create gender map file
echo -e "SAMPLE\tSEX" > "${WORK_DIR}/gender.map"
for ((i=0; i<TOTAL_FILES; i++)); do
    filename=$(basename "${BAM_FILES[$i]}" .bam)
    echo -e "${filename}\tM" >> "${WORK_DIR}/gender.map"  # Default to male
done

# Create ploidy map file
echo -e "SAMPLE\tPLOIDY" > "${WORK_DIR}/ploidy.map"
for ((i=0; i<TOTAL_FILES; i++)); do
    filename=$(basename "${BAM_FILES[$i]}" .bam)
    echo -e "${filename}\t2" >> "${WORK_DIR}/ploidy.map"  # Default to diploid
done

# Create reference chromosome list (exclude unplaced contigs)
echo "Creating reference chromosome list..."
samtools idxstats "${OUTPUT_DIR}/hg38_bams/$(ls -1 ${OUTPUT_DIR}/hg38_bams/ | head -1)" | \
    awk '$1 ~ /^(chr)?[0-9XY]+$/ {print $1}' > "${WORK_DIR}/reference.chrom.list"

# Create metadata file for GenomeSTRIP
echo -e "SAMPLE_ID\tSAMPLE_PATH\tSEQUENCING_CENTER\tPLATFORM\tREFERENCE_GENOME" > "${WORK_DIR}/metadata.txt"
find "${OUTPUT_DIR}/hg38_bams" -name "*.bam" | while read bam_file; do
    filename=$(basename "$bam_file" .hg38.bam)
    echo -e "${filename}\t${bam_file}\tLOCAL\tILLUMINA\thg38" >> "${WORK_DIR}/metadata.txt"
done

# Run SVPreprocess
echo "Running SVPreprocess on all samples..."
cd "$GENOMESTRIP_DIR"
CLASSPATH="${GENOMESTRIP_DIR}/lib/SVToolkit.jar:${GENOMESTRIP_DIR}/lib/gatk/GenomeAnalysisTK.jar:${GENOMESTRIP_DIR}/lib/gatk/Queue.jar"

java -Xmx${MAX_MEMORY} -cp "$CLASSPATH" \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S "${GENOMESTRIP_DIR}/qscript/SVPreprocess.q" \
    -S "${GENOMESTRIP_DIR}/qscript/SVQScript.q" \
    -cp "$CLASSPATH" \
    -gatk "${GENOMESTRIP_DIR}/lib/gatk/GenomeAnalysisTK.jar" \
    -configFile "${GENOMESTRIP_DIR}/conf/genstrip_parameters.txt" \
    -R "$HG38_REFERENCE" \
    -I "${WORK_DIR}/metadata.txt" \
    -md "${OUTPUT_DIR}/genomestrip/preprocess" \
    -bamFilesAreDisjoint true \
    -jobLogDir "${OUTPUT_DIR}/genomestrip/logs" \
    -run

# Run SVDiscovery for deletions
echo "Running SVDiscovery for deletions..."
java -Xmx${MAX_MEMORY} -cp "$CLASSPATH" \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S "${GENOMESTRIP_DIR}/qscript/SVDiscovery.q" \
    -S "${GENOMESTRIP_DIR}/qscript/SVQScript.q" \
    -cp "$CLASSPATH" \
    -gatk "${GENOMESTRIP_DIR}/lib/gatk/GenomeAnalysisTK.jar" \
    -configFile "${GENOMESTRIP_DIR}/conf/genstrip_parameters.txt" \
    -R "$HG38_REFERENCE" \
    -md "${OUTPUT_DIR}/genomestrip/preprocess" \
    -runDirectory "${OUTPUT_DIR}/genomestrip/sv_discovery" \
    -jobLogDir "${OUTPUT_DIR}/genomestrip/logs" \
    -genderMapFile "${WORK_DIR}/gender.map" \
    -ploidyMapFile "${WORK_DIR}/ploidy.map" \
    -I "${WORK_DIR}/metadata.txt" \
    -O "${OUTPUT_DIR}/genomestrip/discovery.deletions.vcf" \
    -minimumSize 100 \
    -maximumSize 1000000 \
    -P deletions \
    -intervalList "${WORK_DIR}/reference.chrom.list" \
    -run

# Run CNVDiscovery for duplications
echo "Running CNVDiscovery pipeline..."
java -Xmx${MAX_MEMORY} -cp "$CLASSPATH" \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S "${GENOMESTRIP_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q" \
    -S "${GENOMESTRIP_DIR}/qscript/SVQScript.q" \
    -cp "$CLASSPATH" \
    -gatk "${GENOMESTRIP_DIR}/lib/gatk/GenomeAnalysisTK.jar" \
    -configFile "${GENOMESTRIP_DIR}/conf/genstrip_parameters.txt" \
    -R "$HG38_REFERENCE" \
    -I "${WORK_DIR}/metadata.txt" \
    -genderMapFile "${WORK_DIR}/gender.map" \
    -md "${OUTPUT_DIR}/genomestrip/preprocess" \
    -runDirectory "${OUTPUT_DIR}/genomestrip/cnv_discovery" \
    -jobLogDir "${OUTPUT_DIR}/genomestrip/logs" \
    -intervalList "${WORK_DIR}/reference.chrom.list" \
    -tilingWindowSize "$WINDOW_SIZE" \
    -tilingWindowOverlap "$WINDOW_OVERLAP" \
    -maximumReferenceGapLength "$REF_GAP_LENGTH" \
    -boundaryPrecision "$BOUNDARY_PRECISION" \
    -minimumRefinedLength "$MIN_REFINED_LENGTH" \
    -O "${OUTPUT_DIR}/genomestrip/discovery.duplications.vcf" \
    -run

echo "STEP 2 completed: GenomeSTRIP analysis"

# ===================== STEP 3: ANNOVAR ANNOTATION =====================
echo "STEP 3: Annotating GenomeSTRIP results with ANNOVAR"

cd "$ANNOVAR_DIR"

# Convert VCF files to ANNOVAR input format
echo "Converting VCF to ANNOVAR format..."
perl convert2annovar.pl -format vcf4 "${OUTPUT_DIR}/genomestrip/discovery.deletions.vcf" > "${OUTPUT_DIR}/annovar/deletions.avinput"
perl convert2annovar.pl -format vcf4 "${OUTPUT_DIR}/genomestrip/discovery.duplications.vcf" > "${OUTPUT_DIR}/annovar/duplications.avinput"

# Create custom annotation database from gEVE GTF file
echo "Creating custom annotation database from gEVE GTF..."
perl gtfToGenePred "$GEVE_GTF" "${OUTPUT_DIR}/annovar/geve_genepred.txt"
perl retrieve_seq_from_fasta.pl -format refGene -seqfile "$HG38_REFERENCE" "${OUTPUT_DIR}/annovar/geve_genepred.txt" -out "${OUTPUT_DIR}/annovar/geve_refGeneMrna.fa"

# Run ANNOVAR annotation with gEVE
echo "Running ANNOVAR annotation..."
perl table_annovar.pl "${OUTPUT_DIR}/annovar/deletions.avinput" \
    "${OUTPUT_DIR}/annovar" \
    -buildver hg38 \
    -out "deletions_annotated" \
    -protocol refGene,geve_genepred \
    -operation g,g \
    -nastring . \
    -remove \
    -otherinfo

perl table_annovar.pl "${OUTPUT_DIR}/annovar/duplications.avinput" \
    "${OUTPUT_DIR}/annovar" \
    -buildver hg38 \
    -out "duplications_annotated" \
    -protocol refGene,geve_genepred \
    -operation g,g \
    -nastring . \
    -remove \
    -otherinfo

echo "STEP 3 completed: ANNOVAR annotation"

# ===================== STEP 4: ERVCALLER ANALYSIS =====================
echo "STEP 4: Running ERVcaller for insertions"

# Set up ERVcaller environment
export PATH=$PATH:${ERVCALLER_DIR}
export PATH=$PATH:${ERVCALLER_DIR}/Scripts

# Create a file listing all BAM files
echo "Creating BAM list for ERVcaller..."
find "${OUTPUT_DIR}/hg38_bams" -name "*.bam" > "${WORK_DIR}/bam_list.txt"

# Run ERVcaller with all BAM files
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
    cp "${OUTPUT_DIR}/ervcaller_analysis/"*".vcf" "${OUTPUT_DIR}/ervcaller_analysis/final-results.vcf"
fi

# Convert VCF to tab-delimited format
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

echo "STEP 4 completed: ERVcaller analysis"

# ===================== STEP 5: INTEGRATION =====================
echo "STEP 5: Integrating and analyzing HERV results"

# Create R script to integrate all data
cat > "${WORK_DIR}/integrate_hervs.R" << 'EOF'
# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(RColorBrewer)
})

# Set working directory and read environment variables
output_dir <- Sys.getenv("OUTPUT_DIR")
geve_gtf <- Sys.getenv("GEVE_GTF")
setwd(output_dir)

# Read gEVE GTF to get HERV reference locations
message("Reading gEVE GTF file...")
geve <- import(geve_gtf)
herv_ref <- geve[geve$type == "gene" | geve$type == "exon"]
herv_ref_gr <- GRanges(
  seqnames = seqnames(herv_ref),
  ranges = ranges(herv_ref),
  strand = strand(herv_ref),
  name = herv_ref$gene_name
)

# Process GenomeSTRIP deletions
message("Processing GenomeSTRIP deletions...")
deletions_file <- "annovar/deletions_annotated.txt"
if (file.exists(deletions_file)) {
  deletions <- read.table(deletions_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="", quote="")
  
  # Create GRanges object
  deletions_gr <- GRanges(
    seqnames = deletions$Chr,
    ranges = IRanges(start = as.numeric(deletions$Start), end = as.numeric(deletions$End)),
    sample = deletions$Otherinfo1,
    gene = deletions$Gene.refGene
  )
  
  # Find overlaps with HERVs
  deletions_overlaps <- findOverlaps(deletions_gr, herv_ref_gr)
  
  # Add HERV annotations
  if (length(deletions_overlaps) > 0) {
    deletions_gr$herv <- NA
    deletions_gr$herv[queryHits(deletions_overlaps)] <- herv_ref_gr$name[subjectHits(deletions_overlaps)]
    
    # Create summary table
    del_summary <- data.frame(
      chromosome = seqnames(deletions_gr),
      start = start(deletions_gr),
      end = end(deletions_gr),
      sample = deletions_gr$sample,
      gene = deletions_gr$gene,
      herv = deletions_gr$herv,
      type = "deletion",
      size = width(deletions_gr)
    )
    
    # Write results
    write.table(del_summary, file="results/herv_deletions.txt", sep="\t", quote=FALSE, row.names=FALSE)
  } else {
    message("No HERV deletions found")
    del_summary <- NULL
  }
} else {
  message("Deletions file not found")
  del_summary <- NULL
}

# Process GenomeSTRIP duplications
message("Processing GenomeSTRIP duplications...")
duplications_file <- "annovar/duplications_annotated.txt"
if (file.exists(duplications_file)) {
  duplications <- read.table(duplications_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="", quote="")
  
  # Create GRanges object
  duplications_gr <- GRanges(
    seqnames = duplications$Chr,
    ranges = IRanges(start = as.numeric(duplications$Start), end = as.numeric(duplications$End)),
    sample = duplications$Otherinfo1,
    gene = duplications$Gene.refGene
  )
  
  # Find overlaps with HERVs
  duplications_overlaps <- findOverlaps(duplications_gr, herv_ref_gr)
  
  # Add HERV annotations
  if (length(duplications_overlaps) > 0) {
    duplications_gr$herv <- NA
    duplications_gr$herv[queryHits(duplications_overlaps)] <- herv_ref_gr$name[subjectHits(duplications_overlaps)]
    
    # Create summary table
    dup_summary <- data.frame(
      chromosome = seqnames(duplications_gr),
      start = start(duplications_gr),
      end = end(duplications_gr),
      sample = duplications_gr$sample,
      gene = duplications_gr$gene,
      herv = duplications_gr$herv,
      type = "duplication",
      size = width(duplications_gr)
    )
    
    # Write results
    write.table(dup_summary, file="results/herv_duplications.txt", sep="\t", quote=FALSE, row.names=FALSE)
  } else {
    message("No HERV duplications found")
    dup_summary <- NULL
  }
} else {
  message("Duplications file not found")
  dup_summary <- NULL
}

# Process ERVcaller insertions
message("Processing ERVcaller insertions...")
insertions_file <- "ervcaller_analysis/final-results.txt"
if (file.exists(insertions_file)) {
  insertions <- read.table(insertions_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  # Create GRanges object
  insertions_gr <- GRanges(
    seqnames = insertions$Chromosome,
    ranges = IRanges(start = as.numeric(insertions$Position), end = as.numeric(insertions$Position) + as.numeric(insertions$Length)),
    sample = insertions$Sample,
    type = insertions$Type,
    herv_type = ifelse(grepl("HERV", insertions$Type), insertions$Type, NA)
  )
  
  # Find overlaps with HERVs
  insertions_overlaps <- findOverlaps(insertions_gr, herv_ref_gr)
  
  # Add HERV annotations
  if (length(insertions_overlaps) > 0) {
    insertions_gr$herv <- NA
    insertions_gr$herv[queryHits(insertions_overlaps)] <- herv_ref_gr$name[subjectHits(insertions_overlaps)]
    
    # Create summary table
    ins_summary <- data.frame(
      chromosome = seqnames(insertions_gr),
      position = start(insertions_gr),
      sample = insertions_gr$sample,
      type = "insertion",
      herv_type = insertions_gr$herv_type,
      herv = insertions_gr$herv,
      size = width(insertions_gr)
    )
    
    # Write results
    write.table(ins_summary, file="results/herv_insertions.txt", sep="\t", quote=FALSE, row.names=FALSE)
  } else {
    message("No HERV insertions found")
    ins_summary <- NULL
  }
} else {
  message("Insertions file not found")
  ins_summary <- NULL
}

# Combine all HERV variants into one dataset
message("Generating combined HERV report...")
all_hervs <- NULL

if (!is.null(del_summary)) {
  all_hervs <- rbind(all_hervs, del_summary[!is.na(del_summary$herv),])
}

if (!is.null(dup_summary)) {
  all_hervs <- rbind(all_hervs, dup_summary[!is.na(dup_summary$herv),])
}

if (!is.null(ins_summary)) {
  # Create compatible data frame structure
  ins_for_merge <- data.frame(
    chromosome = ins_summary$chromosome,
    start = ins_summary$position,
    end = ins_summary$position + ins_summary$size,
    sample = ins_summary$sample,
    gene = NA,
    herv = ins_summary$herv,
    type = ins_summary$type,
    size = ins_summary$size
  )
  all_hervs <- rbind(all_hervs, ins_for_merge[!is.na(ins_for_merge$herv),])
}

if (!is.null(all_hervs) && nrow(all_hervs) > 0) {
  # Write combined results
  write.table(all_hervs, file="results/all_herv_variants.txt", sep="\t", quote=FALSE, row.names=FALSE)
  
  # Generate summary statistics
  herv_counts <- all_hervs %>%
    group_by(herv, type) %>%
    summarize(count = n(), .groups='drop') %>%
    pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric))))
  
  # Write HERV counts
  write.table(herv_counts, file="results/herv_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)
  
  # Generate sample summary
  sample_summary <- all_hervs %>%
    group_by(sample, type) %>%
    summarize(count = n(), .groups='drop') %>%
    pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric))))
  
  # Write sample summary
  write.table(sample_summary, file="results/sample_summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
  
  # Create visualizations
  message("Generating visualizations...")
  
  # Plot HERV counts by type
  if (nrow(herv_counts) > 0) {
    # Prepare data for plotting
    plot_data <- herv_counts %>%
      arrange(desc(total)) %>%
      head(20) %>%
      pivot_longer(cols = c("deletion", "duplication", "insertion"), 
                   names_to = "variant_type", 
                   values_to = "count") %>%
      mutate(herv = factor(herv, levels = herv_counts$herv[order(-herv_counts$total)]))
    
    # Create plot
    herv_plot <- ggplot(plot_data, aes(x = herv, y = count, fill = variant_type)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Top 20 HERVs with Structural Variants",
           x = "HERV",
           y = "Count",
           fill = "Variant Type") +
      scale_fill_brewer(palette = "Set1")
    
    # Save plot
    ggsave("results/top_hervs_plot.pdf", plot = herv_plot, width = 10, height = 6)
  }
  
  # Create comprehensive report
  message("Creating final report...")
  sink("results/herv_analysis_report.txt")
  
  cat("HERV Analysis Report\n")
  cat("===================\n\n")
  
  cat("Summary Statistics:\n")
  cat("-----------------\n")
  cat("Total number of HERV deletions: ", sum(all_hervs$type == "deletion"), "\n")
  cat("Total number of HERV duplications: ", sum(all_hervs$type == "duplication"), "\n")
  cat("Total number of HERV insertions: ", sum(all_hervs$type == "insertion"), "\n\n")
  
  cat("Top HERVs with structural variants:\n")
  print(herv_counts %>% arrange(desc(total)) %>% head(10))
  cat("\n")
  
  cat("Sample-level statistics:\n")
  print(sample_summary)
  cat("\n")
  
  if (sum(all_hervs$type == "deletion") > 0) {
    cat("Deletion size statistics:\n")
    del_sizes <- all_hervs$size[all_hervs$type == "deletion"]
    cat("  Min:", min(del_sizes), "bp\n")
    cat("  Median:", median(del_sizes), "bp\n")
    cat("  Mean:", mean(del_sizes), "bp\n")
    cat("  Max:", max(del_sizes), "bp\n\n")
  }
  
  if (sum(all_hervs$type == "duplication") > 0) {
    cat("Duplication size statistics:\n")
    dup_sizes <- all_hervs$size[all_hervs$type == "duplication"]
    cat("  Min:", min(dup_sizes), "bp\n")
    cat("  Median:", median(dup_sizes), "bp\n")
    cat("  Mean:", mean(dup_sizes), "bp\n")
    cat("  Max:", max(dup_sizes), "bp\n\n")
  }
  
  if (sum(all_hervs$type == "insertion") > 0) {
    cat("Insertion size statistics:\n")
    ins_sizes <- all_hervs$size[all_hervs$type == "insertion"]
    cat("  Min:", min(ins_sizes), "bp\n")
    cat("  Median:", median(ins_sizes), "bp\n")
    cat("  Mean:", mean(ins_sizes), "bp\n")
    cat("  Max:", max(ins_sizes), "bp\n\n")
  }
  
  cat("Chromosomal distribution of HERV variants:\n")
  chrom_dist <- all_hervs %>%
    group_by(chromosome, type) %>%
    summarize(count = n(), .groups='drop')
  print(chrom_dist)
  
  sink()
  
  # Create chromosome distribution plot
  chrom_plot <- ggplot(chrom_dist, aes(x = chromosome, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Chromosomal Distribution of HERV Variants",
         x = "Chromosome",
         y = "Count",
         fill = "Variant Type") +
    scale_fill_brewer(palette = "Set1")
  
  # Save plot
  ggsave("results/chromosome_distribution.pdf", plot = chrom_plot, width = 12, height = 6)
} else {
  message("No HERV variants found in the data")
  sink("results/herv_analysis_report.txt")
  cat("HERV Analysis Report\n")
  cat("===================\n\n")
  cat("No HERV variants were detected in the analyzed data.\n")
  sink()
}

message("Analysis completed. Results are in the 'results' directory.")
EOF

# Run the R integration script
export OUTPUT_DIR="$OUTPUT_DIR"
export GEVE_GTF="$GEVE_GTF"
Rscript "${WORK_DIR}/integrate_hervs.R"

# Generate final summary report
echo "Creating final summary..."

cat > "${OUTPUT_DIR}/results/final_summary.txt" << EOF
======================================================
HERV Analysis Pipeline - Summary Report
======================================================

Date: $(date)

1. BAM Conversion
----------------
Total BAM files processed: $TOTAL_FILES

2. GenomeSTRIP Structural Variant Analysis
----------------------------------------
Deletions: $(grep -c -v "^#" "${OUTPUT_DIR}/genomestrip/discovery.deletions.vcf" 2>/dev/null || echo "0")
Duplications: $(grep -c -v "^#" "${OUTPUT_DIR}/genomestrip/discovery.duplications.vcf" 2>/dev/null || echo "0")

3. ERVcaller Insertion Analysis
----------------------------
Insertions: $(grep -c -v "^#" "${OUTPUT_DIR}/ervcaller_analysis/final-results.vcf" 2>/dev/null || echo "0")

4. HERV Variants Detected
-----------------------
$(if [ -f "${OUTPUT_DIR}/results/herv_counts.txt" ]; then
    echo "Total HERVs affected: $(wc -l < "${OUTPUT_DIR}/results/herv_counts.txt")"
    echo ""
    echo "Top 5 HERVs with variants:"
    head -n 6 "${OUTPUT_DIR}/results/herv_counts.txt"
else
    echo "No HERV variants detected"
fi)

For detailed results, see the analysis report at:
${OUTPUT_DIR}/results/herv_analysis_report.txt
EOF

echo "HERV analysis pipeline completed successfully!"
echo "Final summary is available at ${OUTPUT_DIR}/results/final_summary.txt"
