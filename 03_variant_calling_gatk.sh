#!/bin/bash
#===============================================================================
# STEP 03: Variant Calling - GATK HaplotypeCaller
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="gatk"
log_info "===== STEP 03: ${CALLER} HaplotypeCaller ====="
start_timer

FINAL_BAM="${PREPROC_DIR}/${PREFIX}_final.bam"
check_tool gatk    || exit 1
check_tool bcftools || exit 1
check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"
ensure_dir "${OUT_DIR}"

#-------------------------------------------------------------------------------
# 1. Run HaplotypeCaller
#-------------------------------------------------------------------------------
log_info "Running GATK HaplotypeCaller..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"

run_with_metrics "${CALLER}" "haplotypecaller" "${LOG_DIR}/${CALLER}.log" \
    gatk HaplotypeCaller \
    --java-options "${JAVA_OPTS}" \
    -R "${REF_FASTA}" \
    -I "${FINAL_BAM}" \
    -O "${RAW_VCF}" \
    -L "${WES_BED}" \
    --standard-min-confidence-threshold-for-calling "${GATK_STAND_CALL_CONF}" \
    --native-pair-hmm-threads "${THREADS}"

log_info "HaplotypeCaller completed"

#-------------------------------------------------------------------------------
# 2. Hard filtering (GATK Best Practices)
#-------------------------------------------------------------------------------
log_info "Applying hard filters..."

# SNPs
gatk SelectVariants -V "${RAW_VCF}" -select-type SNP -O "${OUT_DIR}/snps_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/snps_raw.vcf.gz" \
    -O "${OUT_DIR}/snps_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_FILTER"

# INDELs
gatk SelectVariants -V "${RAW_VCF}" -select-type INDEL -O "${OUT_DIR}/indels_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/indels_raw.vcf.gz" \
    -O "${OUT_DIR}/indels_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
    --filter-name "INDEL_FILTER"

# Merge
gatk MergeVcfs \
    -I "${OUT_DIR}/snps_filtered.vcf.gz" \
    -I "${OUT_DIR}/indels_filtered.vcf.gz" \
    -O "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz"

#-------------------------------------------------------------------------------
# 3. Extract PASS variants
#-------------------------------------------------------------------------------
log_info "Extracting PASS variants..."

PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"
bcftools view -f PASS "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz" -Oz -o "${PASS_VCF}"
tabix -p vcf "${PASS_VCF}"

# Split by type for benchmarking
bcftools view -v snps "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${PASS_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Normalize for benchmarking
#-------------------------------------------------------------------------------
log_info "Normalizing variants for benchmarking..."

NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
"${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

#-------------------------------------------------------------------------------
# 5. Stats
#-------------------------------------------------------------------------------
bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${PASS_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${PASS_VCF}" | wc -l)

log_info "Results:  $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

end_timer "03_${CALLER}"
log_info "===== ${CALLER} Complete ====="
