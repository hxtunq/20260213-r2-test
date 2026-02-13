#!/bin/bash
#===============================================================================
# STEP 07: Functional risk assessment of variant-calling errors (FP/FN)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Functional risk assessment ====="
start_timer

check_tool bcftools || exit 1
check_tool tabix || exit 1

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
    ensure_dir "$(dirname "${TRUTH_NORM}")"
    "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

RISK_DIR="${RESULTS_DIR}/functional_risk"
ensure_dir "${RISK_DIR}"
ensure_dir "${RISK_DIR}/errors"

# WES-standard filter (pragmatic clinical shortlist style): QUAL >= 30, DP >= 10 when DP exists.
WES_STD_EXPR='QUAL>=30 && (INFO/DP>=10 || INFO/DP=".")'

for caller in gatk deepvariant strelka2 freebayes; do
    RAW_QUERY_VCF="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.norm.vcf.gz"
    if [[ ! -f "${RAW_QUERY_VCF}" ]]; then
        log_warn "Skip ${caller}: normalized PASS VCF not found (${RAW_QUERY_VCF})"
        continue
    fi

    CALLER_DIR="${VARIANT_DIR}/${caller}"
    WES_STD_QUERY_VCF="${CALLER_DIR}/${PREFIX}_${caller}_wes_standard.pass.norm.vcf.gz"

    # Build WES-standard filtered callset for comparison with raw callset
    bcftools view -i "${WES_STD_EXPR}" -Oz -o "${WES_STD_QUERY_VCF}" "${RAW_QUERY_VCF}"
    tabix -f -p vcf "${WES_STD_QUERY_VCF}"

    for dataset in raw wes_standard; do
        case "${dataset}" in
            raw) QUERY_VCF="${RAW_QUERY_VCF}" ;;
            wes_standard) QUERY_VCF="${WES_STD_QUERY_VCF}" ;;
        esac

        CALLER_ERR_DIR="${RISK_DIR}/errors/${caller}/${dataset}"
        ensure_dir "${CALLER_ERR_DIR}"

        log_info "Extracting FP/FN sets for ${caller} (${dataset})"

        # FN: in truth but missing in caller
        bcftools isec -C -w1 -Oz -o "${CALLER_ERR_DIR}/${PREFIX}_${caller}_${dataset}_fn.vcf.gz" \
            "${TRUTH_NORM}" "${QUERY_VCF}"
        tabix -f -p vcf "${CALLER_ERR_DIR}/${PREFIX}_${caller}_${dataset}_fn.vcf.gz"

        # FP: in caller but not in truth
        bcftools isec -C -w1 -Oz -o "${CALLER_ERR_DIR}/${PREFIX}_${caller}_${dataset}_fp.vcf.gz" \
            "${QUERY_VCF}" "${TRUTH_NORM}"
        tabix -f -p vcf "${CALLER_ERR_DIR}/${PREFIX}_${caller}_${dataset}_fp.vcf.gz"
    done

done

PYTHON_BIN=$(get_python_bin || true)
if [[ -z "${PYTHON_BIN}" ]]; then
    log_warn "Python not available: skipped risk-weighted metric aggregation"
    end_timer "07_functional_risk_assessment"
    exit 0
fi

"${PYTHON_BIN}" "${SCRIPT_DIR}/scripts/risk_weighted_eval.py" \
    --variants-dir "${VARIANT_DIR}" \
    --prefix "${PREFIX}" \
    --errors-dir "${RISK_DIR}/errors" \
    --output-summary "${RISK_DIR}/risk_weighted_summary.tsv" \
    --output-details "${RISK_DIR}/risk_weighted_details.tsv" \
    --alpha-missense "${ALPHAMISSENSE_SCORES}" \
    --alpha-genome "${ALPHAGENOME_SCORES}" \
    --varsage "${VARSAGE_SCORES}"

check_exit "Functional risk metrics"

end_timer "07_functional_risk_assessment"
log_info "===== Functional risk assessment complete ====="
