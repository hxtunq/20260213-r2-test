#!/bin/bash
set -euo pipefail

# cach chay
# chmod +x 09_benchmark_happy.sh
# bash 09_benchmark_happy.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

check_tool docker || exit 1
check_tool bcftools || exit 1

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
ensure_dir "$(dirname "${TRUTH_NORM}")"

if [[ ! -f "${TRUTH_NORM}" ]]; then
  log_info "Normalize truth VCF..."
  "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

for c in "${CALLERS[@]}"; do
  QUERY_NORM="${VARIANT_DIR}/${c}/${PREFIX}_${c}_pass.norm.vcf.gz"
  if [[ ! -f "${QUERY_NORM}" ]]; then
    log_warn "Skip ${c}: missing ${QUERY_NORM}"
    continue
  fi

  OUTDIR="${BENCH_DIR}/happy/${c}"
  rm -rf "${OUTDIR}"
  ensure_dir "${OUTDIR}"
  PREFIX_OUT="${OUTDIR}/${PREFIX}_${c}"

  log_info "hap.py: ${c}"
  docker run --rm \
    --user "$(id -u)":"$(id -g)" \
    --cpus "${THREADS}" \
    --memory "${MAX_MEMORY}" \
    -v "${PROJECT_DIR}:/work" \
    -w /work \
    "${HAP_IMAGE}" \
    hap.py \
      "${TRUTH_NORM}" \
      "${QUERY_NORM}" \
      --reference "${REF_FASTA}" \
      --threads "${THREADS}" \
      -f "${HIGH_CONF_BED}" \
      --engine vcfeval \
      -o "${PREFIX_OUT}"

  log_info "Done: ${OUTDIR} (check *.summary.csv)"
done

log_info "All done (hap.py)."
