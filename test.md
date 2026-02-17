```bash
cd /workspace/20260213-r2-test
source config/config.sh
CALLER=gatk
OUT_DIR="${VARIANT_DIR}/${CALLER}"
PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"
NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"

bash scripts/normalize_vcf.sh "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

TRUTH_NORM="${BENCH_DIR}/truth/${PREFIX}_truth.norm.vcf.gz"
if [[ ! -f "${TRUTH_NORM}" ]]; then
  mkdir -p "$(dirname "${TRUTH_NORM}")"
  bash scripts/normalize_vcf.sh "${TRUTH_VCF}" "${TRUTH_NORM}" "${REF_FASTA}"
fi

bcftools stats "${PASS_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"
```
