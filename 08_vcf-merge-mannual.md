```bash
OUTDIR=results/benchmarks/rtg_vcfeval/merged
mkdir -p "$OUTDIR"

FN_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fn.vcf.gz | head -n 1)
FN_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fn.vcf.gz | head -n 1)
FN_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fn.vcf.gz | head -n 1)
FN_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fn.vcf.gz | head -n 1)

bcftools reheader -s <(echo gatk) "$FN_GATK" \
  | bcftools +setGT - -Oz -o "$OUTDIR/gatk.fn.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo deepvariant) "$FN_DV" \
  | bcftools +setGT - -Oz -o "$OUTDIR/deepvariant.fn.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo strelka2) "$FN_ST" \
  | bcftools +setGT - -Oz -o "$OUTDIR/strelka2.fn.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo freebayes) "$FN_FB" \
  | bcftools +setGT - -Oz -o "$OUTDIR/freebayes.fn.vcf.gz" -- -t a -n 1/1

bcftools merge -m none -Oz -o "$OUTDIR/fn_4callers.vcf.gz" \
  "$OUTDIR/gatk.fn.vcf.gz" "$OUTDIR/deepvariant.fn.vcf.gz" \
  "$OUTDIR/strelka2.fn.vcf.gz" "$OUTDIR/freebayes.fn.vcf.gz"

tabix -f -p vcf "$OUTDIR/fn_4callers.vcf.gz"
```

```bash
FP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fp.vcf.gz | head -n 1)
FP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fp.vcf.gz | head -n 1)
FP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fp.vcf.gz | head -n 1)
FP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fp.vcf.gz | head -n 1)

bcftools reheader -s <(echo gatk) "$FP_GATK" \
  | bcftools +setGT - -Oz -o "$OUTDIR/gatk.fp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo deepvariant) "$FP_DV" \
  | bcftools +setGT - -Oz -o "$OUTDIR/deepvariant.fp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo strelka2) "$FP_ST" \
  | bcftools +setGT - -Oz -o "$OUTDIR/strelka2.fp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo freebayes) "$FP_FB" \
  | bcftools +setGT - -Oz -o "$OUTDIR/freebayes.fp.vcf.gz" -- -t a -n 1/1

bcftools merge -m none -Oz -o "$OUTDIR/fp_4callers.vcf.gz" \
  "$OUTDIR/gatk.fp.vcf.gz" "$OUTDIR/deepvariant.fp.vcf.gz" \
  "$OUTDIR/strelka2.fp.vcf.gz" "$OUTDIR/freebayes.fp.vcf.gz"

tabix -f -p vcf "$OUTDIR/fp_4callers.vcf.gz"
```

```bash
TP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*tp.vcf.gz | head -n 1)
TP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*tp.vcf.gz | head -n 1)
TP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*tp.vcf.gz | head -n 1)
TP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*tp.vcf.gz | head -n 1)

bcftools reheader -s <(echo gatk) "$TP_GATK" \
  | bcftools +setGT - -Oz -o "$OUTDIR/gatk.tp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo deepvariant) "$TP_DV" \
  | bcftools +setGT - -Oz -o "$OUTDIR/deepvariant.tp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo strelka2) "$TP_ST" \
  | bcftools +setGT - -Oz -o "$OUTDIR/strelka2.tp.vcf.gz" -- -t a -n 1/1
bcftools reheader -s <(echo freebayes) "$TP_FB" \
  | bcftools +setGT - -Oz -o "$OUTDIR/freebayes.tp.vcf.gz" -- -t a -n 1/1

bcftools merge -m none -Oz -o "$OUTDIR/tp_4callers.vcf.gz" \
  "$OUTDIR/gatk.tp.vcf.gz" "$OUTDIR/deepvariant.tp.vcf.gz" \
  "$OUTDIR/strelka2.tp.vcf.gz" "$OUTDIR/freebayes.tp.vcf.gz"

tabix -f -p vcf "$OUTDIR/tp_4callers.vcf.gz"
```
