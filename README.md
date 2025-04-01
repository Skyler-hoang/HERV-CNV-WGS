python ervcaller_to_jasmine.py -i input.vcf -o output_dir -g 20


# HERV-CNV-WGS
# HPC Usage Guide for HERV Analysis Pipeline

This guide explains how to use the HERV analysis pipeline on an HPC system.

## Prerequisites

Before running the script, ensure the following tools are installed on your HPC:

- BWA (version 0.7.10 or later)
- Samtools (version 1.6 or later)
- CrossMap
- GenomeSTRIP
- ANNOVAR
- ERVcaller (version 1.4)
- R (with GenomicRanges, rtracklayer, ggplot2, dplyr, tidyr packages)

## Setting Up the Environment

Many HPC systems use module systems to load software. Check if required tools are available:

```bash
module avail
```

Load required modules:

```bash
module load bwa
module load samtools
module load r
# Load any other available modules
```

For tools not available as modules, set up your environment:

```bash
# Example for ERVcaller
export PATH=$PATH:/path/to/ERVcaller_v1.4
```

## Step 1: Prepare the Script

1. Download the script:
   ```bash
   wget https://your-script-url/herv_analysis_pipeline.sh
   ```

2. Make the script executable:
   ```bash
   chmod +x herv_analysis_pipeline.sh
   ```

3. Edit the configuration section at the top of the script:
   ```bash
   nano herv_analysis_pipeline.sh
   ```

   Update these paths to match your HPC environment:
   
   ```bash
   WORK_DIR="/path/to/your/working/directory"
   INPUT_DIR="/path/to/bam/files"
   HG19_REFERENCE="/path/to/references/hg19.fa"
   HG38_REFERENCE="/path/to/references/GRCh38.fa"
   CHAIN_FILE="/path/to/references/hg19ToHg38.over.chain"
   GENOMESTRIP_DIR="/path/to/genomestrip"
   ANNOVAR_DIR="/path/to/annovar"
   ERVCALLER_DIR="/path/to/ERVcaller_v1.4"
   GEVE_GTF="/path/to/gEVE.gtf"
   GEVE_FA="/path/to/gEVE.fa"
   ```

## Step 2: Create a Job Submission Script

Create a file named `submit_job.sh` with the appropriate settings for your HPC system:

### For SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=HERV_analysis
#SBATCH --output=HERV_analysis_%j.out
#SBATCH --error=HERV_analysis_%j.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Load modules
module load bwa
module load samtools
module load r

# Run the analysis script
./herv_analysis_pipeline.sh
```

### For PBS/Torque:

```bash
#!/bin/bash
#PBS -N HERV_analysis
#PBS -o HERV_analysis.out
#PBS -e HERV_analysis.err
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb

# Load modules
module load bwa
module load samtools
module load r

# Change to submission directory
cd $PBS_O_WORKDIR

# Run the analysis script
./herv_analysis_pipeline.sh
```

## Step 3: Submit the Job

1. Make the submission script executable:
   ```bash
   chmod +x submit_job.sh
   ```

2. Submit the job:

   - For SLURM:
     ```bash
     sbatch submit_job.sh
     ```

   - For PBS/Torque:
     ```bash
     qsub submit_job.sh
     ```

## Step 4: Monitor Job Progress

1. Check job status:
   - SLURM: `squeue -u your_username`
   - PBS: `qstat -u your_username`

2. Monitor log output:
   ```bash
   tail -f herv_analysis.log
   ```

3. Check resource usage:
   - SLURM: `sstat --format=AveCPU,AveRSS,AveVMSize,MaxRSS,MaxVMSize -j job_id`
   - PBS: `qstat -f job_id`

## Step 5: Troubleshooting Common Issues

1. **Out of memory errors**: Increase memory allocation in job script.

2. **GenomeSTRIP failures**: 
   - Ensure Java version is compatible (Java 8 recommended)
   - Check that reference genome is properly indexed
   - Increase Java heap size by modifying MAX_MEMORY variable

3. **ERVcaller issues**:
   - Verify that the BAM list file has correct paths
   - Check that all required dependencies are installed

4. **R package errors**:
   - Install missing R packages:
     ```R
     install.packages(c("GenomicRanges", "rtracklayer", "ggplot2", "dplyr", "tidyr", "gridExtra", "RColorBrewer"))
     ```

## Step 6: Accessing Results

After job completion, results will be in the specified output directory:

```bash
cd /path/to/your/working/directory/output/results
```

Main result files:
- `final_summary.txt`: Overview of detected variants
- `herv_analysis_report.txt`: Detailed analysis results
- `all_herv_variants.txt`: Combined list of all HERV variants
- Various PDF plots showing distributions of detected variants

## Additional Notes for gEVE Files

The gEVE GTF and FASTA files can be prepared as follows:

1. gEVE GTF file format should look like:
   ```
   chr1  gEVE  gene  74855  76774  .  +  .  gene_id "Hsap38.chr1.74855.76774.+"; gene_name "L1PA2|LINE/L1";
   ```

2. If you need to convert the gEVE format to GTF:
   ```bash
   awk '{print $2"\tgEVE\tgene\t"$3"\t"$4"\t.\t"$5"\t.\tgene_id \""$1"\"; gene_name \""$NF"\";";}' gEVE_file.txt > gEVE.gtf
   ```

If you have questions or encounter issues, contact your HPC administrator for assistance.
