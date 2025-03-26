#!/bin/bash
#SBATCH --cpus-per-task=16      # Increase CPU cores (16 instead of 1)
#SBATCH --mem=128G               # Increase memory (64GB instead of 1GB)
#SBATCH --time=24:00:00         # Increase time limit if need

export PATH=$PATH:/home/hhh38/herv_analysis/tools/ERVcaller/Scripts/SE-MEI
export PATH=$PATH:/home/hhh38/bwa

module load R
  
samtools index /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder/CB80_P00050018A.hg38.bam /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder/CB80_P00050018A.hg38.bam.bai

perl /home/hhh38/herv_analysis/tools/ERVcaller/ERVcaller_v1.4.pl \
  -i  CB80_P00050018A.hg38 \
  -f .bam \
  -H /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/reference/hg38/hg38.fa \
  -T /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/reference/TE_consensus.fa \
  -I /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/hg19_bams/temp_placeholder \
  -O /mnt/pan/SOM_EPBI_SKI/PHONOLOGY/Retrotransposons/retro_analysis/output/ervcaller_analysis/ \
  -t 8 \
  -S 20 \
  -G
