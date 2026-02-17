#!/bin/bash
#===============================================================================
# STEP 07: Benchmark PASS+normalized VCF của nhiều caller bằng RTG vcfeval
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/../scripts/helper_functions.sh"

log_info "Start benchmark RTG vcfeval"

missing_tools=()
for tool in docker bcftools; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    missing_tools+=("${tool}")
  fi
done

if [[ ${#missing_tools[@]} -gt 0 ]]; then
  missing_csv=$(IFS=, ; echo "${missing_tools[*]}")
  # In some wrappers stderr is hidden; print to both streams for easier troubleshooting.
  log_error "Missing required tool(s): ${missing_csv}"
  echo "[ERROR] Missing required tool(s): ${missing_csv}"
  exit 1
fi

# Có thể override bằng biến môi trường, ví dụ:
# CALLERS="gatk deepvariant" bash stages/07_benchmark_rtg_vcfeval.sh
CALLERS_RAW="${CALLERS:-gatk deepvariant strelka2 freebayes}"
IFS=' ' read -r -a CALLER_LIST <<< "${CALLERS_RAW}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
SDF="${BENCH_DIR}/ref/${CHR_TO_USE}.sdf"
OUT_ROOT="${BENCH_DIR}/rtg_vcfeval"
SUMMARY_TSV="${OUT_ROOT}/vcfeval_summary.tsv"

ensure_dir "${BENCH_DIR}"
ensure_dir "$(dirname "${TRUTH_NORM}")"
ensure_dir "$(dirname "${SDF}")"
ensure_dir "${OUT_ROOT}"

if [[ ! -f "${TRUTH_VCF}" ]]; then
  log_error "Missing truth VCF: ${TRUTH_VCF}"
  exit 1
fi

if [[ ! -f "${HIGH_CONF_BED}" ]]; then
  log_error "Missing high-confidence BED: ${HIGH_CONF_BED}"
  exit 1
fi

# 1) Normalize truth nếu chưa có
if [[ ! -f "${TRUTH_NORM}" ]]; then
  log_info "Normalize truth VCF -> ${TRUTH_NORM}"
  "${SCRIPT_DIR}/../scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

# 2) Tạo RTG template SDF một lần
if [[ ! -d "${SDF}" ]]; then
  log_info "Create RTG SDF template: ${SDF}"
  docker run --rm \
    --user "$(id -u)":"$(id -g)" \
    -v "${PROJECT_DIR}:/work" \
    -w /work \
    "${RTG_IMAGE}" \
    rtg format --output "${SDF}" "${REF_FASTA}"
fi

# Header summary
printf "caller\ttp\tfp\tfn\tprecision\trecall\tf1\n" > "${SUMMARY_TSV}"

# 3) Chạy vcfeval cho từng caller
for caller in "${CALLER_LIST[@]}"; do
  QUERY_NORM="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.norm.vcf.gz"

  if [[ ! -f "${QUERY_NORM}" ]]; then
    log_warn "Skip ${caller}: missing ${QUERY_NORM}"
    continue
  fi

  OUTDIR="${OUT_ROOT}/${caller}"
  rm -rf "${OUTDIR}"

  log_info "Run RTG vcfeval for caller=${caller}"
  docker run --rm \
    --user "$(id -u)":"$(id -g)" \
    --cpus "${THREADS}" \
    --memory "${MAX_MEMORY}" \
    -v "${PROJECT_DIR}:/work" \
    -w /work \
    "${RTG_IMAGE}" \
    env RTG_MEM="${MAX_MEMORY}" rtg vcfeval \
      --baseline "${TRUTH_NORM}" \
      --calls "${QUERY_NORM}" \
      --template "${SDF}" \
      --bed-regions "${HIGH_CONF_BED}" \
      --output "${OUTDIR}" \
      --threads "${THREADS}"

  SUMMARY_FILE="${OUTDIR}/summary.txt"
  if [[ -f "${SUMMARY_FILE}" ]]; then
    # Parse robust theo format summary.txt của RTG
    tp=$(awk -F':' '/True Positives/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")
    fp=$(awk -F':' '/False Positives/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")
    fn=$(awk -F':' '/False Negatives/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")
    precision=$(awk -F':' '/Precision/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")
    recall=$(awk -F':' '/Sensitivity/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")
    f1=$(awk -F':' '/F-measure/{gsub(/ /, "", $2); print $2; exit}' "${SUMMARY_FILE}" || echo "NA")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${caller}" "${tp}" "${fp}" "${fn}" "${precision}" "${recall}" "${f1}" >> "${SUMMARY_TSV}"
  else
    log_warn "Missing summary.txt for ${caller}: ${SUMMARY_FILE}"
  fi

done

log_info "Done. Combined summary: ${SUMMARY_TSV}"
