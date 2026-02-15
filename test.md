```bash
REF=data/reference/chr22.fa
BED=data/reference/chr22_non_N_regions.bed
SDF=results/benchmarks/ref/chr22.sdf
BASELINE=results/benchmarks/truth/truth.gt.norm.vcf.gz

GATK_FILT=$(ls -1 results/variants/gatk/*_gatk_filtered.vcf.gz | head -n 1)

# 1) Lấy ALL (không -f)
ALL_VCF=results/variants/gatk/gatk_all.vcf.gz
bcftools view "$GATK_FILT" -Oz -o "$ALL_VCF"
tabix -f -p vcf "$ALL_VCF"

# 2) Normalize ALL
ALL_NORM=results/variants/gatk/gatk_all.norm.vcf.gz
TMP=${ALL_NORM}.tmp.vcf.gz
bcftools norm -f "$REF" -m -both "$ALL_VCF" -Oz -o "$TMP"
bcftools sort "$TMP" -Oz -o "$ALL_NORM"
rm -f "$TMP"
tabix -f -p vcf "$ALL_NORM"

OUTDIR=results/benchmarks/rtg_vcfeval/gatk_all

rm -rf "$OUTDIR"
mkdir -p "$(dirname "$OUTDIR")"

RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls results/variants/gatk/gatk_all.norm.vcf.gz \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions data/reference/chr22_non_N_regions.bed \
  --output "$OUTDIR" \
  --threads 4

cat "$OUTDIR"/*.summary.txt

```
