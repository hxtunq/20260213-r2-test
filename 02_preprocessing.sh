#!/bin/bash
#===============================================================================
# STEP 02: Preprocessing (GATK Best Practices - nf-core/sarek)
# Pipeline: FastQC → BWA-MEM → MarkDuplicates → BQSR
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 02: Preprocessing (GATK Best Practices) ====="
start_timer

# Input
R1="${SIM_DIR}/${PREFIX}_R1.fastq.gz"
R2="${SIM_DIR}/${PREFIX}_R2.fastq.gz"

check_file "${R1}" || exit 1
check_file "${R2}" || exit 1

#-------------------------------------------------------------------------------
# 1. FastQC - Raw reads quality control
#-------------------------------------------------------------------------------
log_info "[1/4] FastQC on raw reads..."

ensure_dir "${PREPROC_DIR}/fastqc_raw"
fastqc -t "${THREADS}" -o "${PREPROC_DIR}/fastqc_raw" "${R1}" "${R2}" \
    2>&1 | tee "${LOG_DIR}/fastqc_raw.log"

check_exit "FastQC raw"

#-------------------------------------------------------------------------------
# 2. BWA-MEM - Alignment (directly from raw FASTQ after FastQC)
#-------------------------------------------------------------------------------
log_info "[2/4] BWA-MEM alignment..."

ALIGNED_BAM="${PREPROC_DIR}/${PREFIX}_aligned.bam"

bwa mem \
    -t "${THREADS}" \
    -R "${READ_GROUP}" \
    -M \
    "${REF_FASTA}" \
    "${R1}" "${R2}" \
    2> "${LOG_DIR}/bwa_mem.log" | \
samtools sort -@ "${THREADS}" -m 2G -o "${ALIGNED_BAM}" -

samtools index "${ALIGNED_BAM}"

check_exit "BWA-MEM"

#-------------------------------------------------------------------------------
# 3. GATK MarkDuplicates
#-------------------------------------------------------------------------------
log_info "[3/4] GATK MarkDuplicates..."

MARKED_BAM="${PREPROC_DIR}/${PREFIX}_marked.bam"
DUP_METRICS="${PREPROC_DIR}/${PREFIX}_dup_metrics.txt"

gatk MarkDuplicates \
    --java-options "${JAVA_OPTS}" \
    -I "${ALIGNED_BAM}" \
    -O "${MARKED_BAM}" \
    -M "${DUP_METRICS}" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --CREATE_INDEX true \
    2>&1 | tee "${LOG_DIR}/markduplicates.log"

check_exit "MarkDuplicates"

#-------------------------------------------------------------------------------
# 4. GATK BaseRecalibrator & ApplyBQSR
#-------------------------------------------------------------------------------
log_info "[4/4] Base Quality Score Recalibration..."

# Ensure chromosome naming consistency
ensure_chr_naming() {
    local vcf_file="$1"
    if [[ ! -f "${vcf_file}" ]]; then
        return 1
    fi
    
    local first_chr=$(bcftools view -H "${vcf_file}" 2>/dev/null | head -1 | cut -f1)
    if [[ -n "${first_chr}" && "${first_chr}" != chr* ]]; then
        log_warn "Renaming chromosomes in ${vcf_file} to match reference..."
        local tmp_renamed="${vcf_file%.vcf.gz}_renamed.vcf.gz"
        "${SCRIPT_DIR}/scripts/rename_chromosomes.sh" "${vcf_file}" "${tmp_renamed}"
        mv "${tmp_renamed}" "${vcf_file}"
        mv "${tmp_renamed}.tbi" "${vcf_file}.tbi" 2>/dev/null || tabix -p vcf "${vcf_file}"
    fi
    return 0
}

RECAL_TABLE="${PREPROC_DIR}/${PREFIX}_recal.table"
FINAL_BAM="${PREPROC_DIR}/${PREFIX}_recal.bam"

log_info "Running BQSR with known sites..."

# Ensure chromosome naming convention matches reference (chr22)
# This is critical for BQSR to work correctly
for vcf in "${DBSNP}" "${KNOWN_INDELS}"; do
    if [[ -f "${vcf}" ]]; then
        ensure_chr_naming "${vcf}" || log_warn "Failed to verify chromosome naming for ${vcf}"
    fi
done

# BaseRecalibrator
KNOWN_SITES_ARGS=""
[[ -n "${DBSNP:-}" ]] && KNOWN_SITES_ARGS+="--known-sites ${DBSNP} "
[[ -n "${KNOWN_INDELS:-}" ]] && KNOWN_SITES_ARGS+="--known-sites ${KNOWN_INDELS} "

gatk BaseRecalibrator \
    --java-options "${JAVA_OPTS}" \
    -R "${REF_FASTA}" \
    -I "${MARKED_BAM}" \
    ${KNOWN_SITES_ARGS} \
    -L "${WES_BED}" \
    -O "${RECAL_TABLE}" \
    2>&1 | tee "${LOG_DIR}/baserecalibrator.log"

# ApplyBQSR
gatk ApplyBQSR \
    --java-options "${JAVA_OPTS}" \
    -R "${REF_FASTA}" \
    -I "${MARKED_BAM}" \
    --bqsr-recal-file "${RECAL_TABLE}" \
    -L "${WES_BED}" \
    -O "${FINAL_BAM}" \
    2>&1 | tee "${LOG_DIR}/applybqsr.log"

samtools index "${FINAL_BAM}"
check_exit "BQSR"

#-------------------------------------------------------------------------------
# 5. Filter by mapping quality
#-------------------------------------------------------------------------------
log_info "Filtering by mapping quality (MAPQ >= ${MIN_MAPPING_QUALITY})..."

FILTERED_BAM="${PREPROC_DIR}/${PREFIX}_final.bam"

samtools view \
    -@ "${THREADS}" \
    -b \
    -q "${MIN_MAPPING_QUALITY}" \
    -F 1796 \
    "${FINAL_BAM}" | \
samtools sort -@ "${THREADS}" -o "${FILTERED_BAM}" -

samtools index "${FILTERED_BAM}"

check_exit "BAM filtering"

#-------------------------------------------------------------------------------
# 6. Alignment statistics (samtools stats, mosdepth)
#-------------------------------------------------------------------------------
log_info "Calculating alignment statistics..."

samtools stats "${FILTERED_BAM}" > "${PREPROC_DIR}/${PREFIX}_stats.txt"
samtools flagstat "${FILTERED_BAM}" > "${PREPROC_DIR}/${PREFIX}_flagstat.txt"
samtools idxstats "${FILTERED_BAM}" > "${PREPROC_DIR}/${PREFIX}_idxstats.txt"

# Mosdepth for coverage
if check_tool mosdepth; then
    mosdepth -t "${THREADS}" --by 1000 \
        "${PREPROC_DIR}/${PREFIX}_coverage" "${FILTERED_BAM}"
fi

#-------------------------------------------------------------------------------
# 7. Summary
#-------------------------------------------------------------------------------
MAPPED=$(grep "reads mapped:" "${PREPROC_DIR}/${PREFIX}_stats.txt" | cut -f3)
DUP_RATE=$(grep "PERCENT_DUPLICATION" "${DUP_METRICS}" -A1 | tail -1 | cut -f9 || echo "N/A")

log_info "===== Preprocessing Summary ====="
log_info "  Final BAM: ${FILTERED_BAM}"
log_info "  Mapped reads: ${MAPPED}"
log_info "  Duplication rate: ${DUP_RATE}"

# Export for downstream scripts
export FINAL_BAM="${FILTERED_BAM}"
echo "FINAL_BAM=${FILTERED_BAM}" > "${PREPROC_DIR}/bam_path.sh"

end_timer "02_preprocessing"
log_info "===== Preprocessing Complete ====="
