#!/bin/bash
#===============================================================================
# STEP 02: Preprocessing (GATK Best Practices - nf-core/sarek)
# Manual commands without variable expansion.
# Pipeline: FastQC → BWA-MEM → MarkDuplicates → BQSR
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# 1. FastQC - Raw reads quality control
#-------------------------------------------------------------------------------
mkdir -p results/preprocessing/fastqc_raw
fastqc -t 4 \
  -o results/preprocessing/fastqc_raw \
  data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  2>&1 | tee logs/fastqc_raw.log

#-------------------------------------------------------------------------------
# 2. BWA-MEM - Alignment (directly from raw FASTQ after FastQC)
#-------------------------------------------------------------------------------
bwa mem \
  -t 4 \
  -R "@RG\tID:SIMULATED_SAMPLE\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  -M \
  data/reference/chr22.fa \
  data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  2> logs/bwa_mem.log | \
  samtools sort -@ 4 -m 2G -o results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam

#-------------------------------------------------------------------------------
# 3. GATK MarkDuplicates
#-------------------------------------------------------------------------------
gatk MarkDuplicates \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  -M results/preprocessing/SIMULATED_SAMPLE_chr22_dup_metrics.txt \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX true \
  2>&1 | tee logs/markduplicates.log

#-------------------------------------------------------------------------------
# 4. GATK BaseRecalibrator & ApplyBQSR
#-------------------------------------------------------------------------------
# Ensure known-sites VCFs match the "chr" naming convention if needed.
# If your VCFs use "22" instead of "chr22", run:
# scripts/rename_chromosomes.sh <input.vcf.gz> <output.vcf.gz>

gatk BaseRecalibrator \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --known-sites data/reference/dbsnp138.hg38.chr22.vcf.gz \
  --known-sites data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  --known-sites data/reference/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  2>&1 | tee logs/baserecalibrator.log

gatk ApplyBQSR \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --bqsr-recal-file results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam \
  2>&1 | tee logs/applybqsr.log

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam

#-------------------------------------------------------------------------------
# 5. Filter by mapping quality
#-------------------------------------------------------------------------------
samtools view \
  -@ 4 \
  -b \
  -q 20 \
  -F 1796 \
  results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam | \
  samtools sort -@ 4 -o results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam

#-------------------------------------------------------------------------------
# 6. Alignment statistics (samtools stats, mosdepth)
#-------------------------------------------------------------------------------
samtools stats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_stats.txt
samtools flagstat results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_flagstat.txt
samtools idxstats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_idxstats.txt

# mosdepth for coverage
mosdepth -t 4 --by 1000 \
  results/preprocessing/SIMULATED_SAMPLE_chr22_coverage \
  results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam

#-------------------------------------------------------------------------------
# 7. Export for downstream scripts
#-------------------------------------------------------------------------------
echo "FINAL_BAM=results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam" > results/preprocessing/bam_path.sh
