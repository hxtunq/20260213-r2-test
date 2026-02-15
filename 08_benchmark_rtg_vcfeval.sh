#!/bin/bash
set -euo pipefail

# cach chay:
# chmod +x 08_benchmark_rtg_vcfeval.sh
# bash 08_benchmark_rtg_vcfeval.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

check_tool docker || exit 1
check_tool bcftools || exit 1

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
SDF_DIR="${BENCH_DIR}/ref"
SDF="${SDF_DIR}/${CHR_TO_USE}.sdf"

ensure_dir "${BENCH_DIR}"
ensure_dir "$(dirname "${TRUTH_NORM}")"
ensure_dir "${SDF_DIR}"

# Ensure truth normalized (repo bạn đã normalize trong các bước caller, nhưng phòng khi bạn chạy riêng)
if [[ ! -f "${TRUTH_NORM}" ]]; then
  log_info "Normalize truth VCF..."
  "${SCRIPT_DIR}/scripts/normalize_vcf.sh" "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

# Build SDF once
if [[ ! -d "${SDF}" ]]; then
  log_info "Create RTG SDF template: ${SDF}"
  docker run --rm \
    --user "$(id -u)":"$(id -g)" \
    -v "${PROJECT_DIR}:/work" \
    -w /work \
    "${RTG_IMAGE}" \
    rtg format --output "${SDF}" "${REF_FASTA}"
fi

for c in "${CALLERS[@]}"; do
  QUERY_NORM="${VARIANT_DIR}/${c}/${PREFIX}_${c}_pass.norm.vcf.gz"
  if [[ ! -f "${QUERY_NORM}" ]]; then
    log_warn "Skip ${c}: missing ${QUERY_NORM}"
    continue
  fi

  OUTDIR="${BENCH_DIR}/rtg_vcfeval/${c}"
  rm -rf "${OUTDIR}"
  ensure_dir "$(dirname "${OUTDIR}")"

  log_info "RTG vcfeval: ${c}"
  docker run --rm \
    --user "$(id -u)":"$(id -g)" \
    --cpus "${THREADS}" \
    --memory "${MAX_MEMORY}" \
    -e RTG_MEM="${MAX_MEMORY}" \
    -v "${PROJECT_DIR}:/work" \
    -w /work \
    "${RTG_IMAGE}" \
    rtg vcfeval \
      --baseline "${TRUTH_NORM}" \
      --calls "${QUERY_NORM}" \
      --template "${SDF}" \
      --bed-regions "${HIGH_CONF_BED}" \
      --output "${OUTDIR}" \
      --threads "${THREADS}"

  log_info "Done: ${OUTDIR} (check *.summary.txt)"
done

log_info "All done (rtg vcfeval)."
