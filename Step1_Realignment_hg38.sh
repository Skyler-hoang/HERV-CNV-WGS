#!/bin/bash
#SBATCH --cpus-per-task=32      # Increase CPU cores (16 instead of 1)
#SBATCH --mem=128G               # Increase memory (64GB instead of 1GB)
#SBATCH --time=48:00:00         # Increase time limit if need


samtools collate -Oun128 /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder/CB80_P00050018A.bam | samtools fastq -OT RG,BC - \
  | /home/hhh38/bwa/bwa mem -pt8 -CH <(samtools view -H /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder/CB80_P00050018A.bam|grep ^@RG) /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/reference/hg38/hg38.fa - \
  | samtools sort -@4 -m4g -o /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder/CB80_P00050018A.hg38.bam -
