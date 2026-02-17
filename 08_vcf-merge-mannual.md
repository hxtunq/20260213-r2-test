Hợp nhất các file output fn/fp/tp sau khi dùng rtgtool/vcfeval của 4 công cụ gọi biến thể:

``` bash
OUTDIR=results/benchmarks/rtg_vcfeval/merged
mkdir -p "$OUTDIR"
```

```bash
# FN

FN_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fn*.vcf.gz | head -n 1)
FN_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fn*.vcf.gz | head -n 1)
FN_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fn*.vcf.gz | head -n 1)
FN_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fn*.vcf.gz | head -n 1)

for X in gatk:"$FN_GATK" deepvariant:"$FN_DV" strelka2:"$FN_ST" freebayes:"$FN_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.fn.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.fn.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/fn_4callers.vcf.gz" \
  "$OUTDIR/gatk.fn.vcf.gz" "$OUTDIR/deepvariant.fn.vcf.gz" \
  "$OUTDIR/strelka2.fn.vcf.gz" "$OUTDIR/freebayes.fn.vcf.gz"
tabix -f -p vcf "$OUTDIR/fn_4callers.vcf.gz"
```

```bash
# FP

FP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*fp*.vcf.gz | head -n 1)
FP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*fp*.vcf.gz | head -n 1)
FP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*fp*.vcf.gz | head -n 1)
FP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*fp*.vcf.gz | head -n 1)

for X in gatk:"$FP_GATK" deepvariant:"$FP_DV" strelka2:"$FP_ST" freebayes:"$FP_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.fp.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.fp.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/fp_4callers.vcf.gz" \
  "$OUTDIR/gatk.fp.vcf.gz" "$OUTDIR/deepvariant.fp.vcf.gz" \
  "$OUTDIR/strelka2.fp.vcf.gz" "$OUTDIR/freebayes.fp.vcf.gz"
tabix -f -p vcf "$OUTDIR/fp_4callers.vcf.gz"
```

```bash
# TP (TP này là TP_callset, không phải baseline)

TP_GATK=$(ls results/benchmarks/rtg_vcfeval/gatk/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_DV=$(ls results/benchmarks/rtg_vcfeval/deepvariant/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_ST=$(ls results/benchmarks/rtg_vcfeval/strelka2/*tp*.vcf.gz | grep -v baseline | head -n 1)
TP_FB=$(ls results/benchmarks/rtg_vcfeval/freebayes/*tp*.vcf.gz | grep -v baseline | head -n 1)

for X in gatk:"$TP_GATK" deepvariant:"$TP_DV" strelka2:"$TP_ST" freebayes:"$TP_FB"; do
  S=${X%%:*}; IN=${X#*:}
  zcat "$IN" | awk -v s="$S" 'BEGIN{OFS="\t"; gt=0}
    /^##FORMAT=<ID=GT/ {gt=1}
    /^##/ {print; next}
    /^#CHROM/ {
      if(!gt) print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",s
      next
    }
    {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","1/1"}
  ' | bgzip -c > "$OUTDIR/$S.tp.vcf.gz"
  tabix -f -p vcf "$OUTDIR/$S.tp.vcf.gz"
done

bcftools merge -m none -Oz -o "$OUTDIR/tp_4callers.vcf.gz" \
  "$OUTDIR/gatk.tp.vcf.gz" "$OUTDIR/deepvariant.tp.vcf.gz" \
  "$OUTDIR/strelka2.tp.vcf.gz" "$OUTDIR/freebayes.tp.vcf.gz"
tabix -f -p vcf "$OUTDIR/tp_4callers.vcf.gz"
```

Chuyển file VCF sang CSV để xử lý dữ liệu.

```bash
# FN

IN=results/benchmarks/rtg_vcfeval/merged/fn_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/fn_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```

```bash
# FP

IN=results/benchmarks/rtg_vcfeval/merged/fp_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/fp_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```

```bash
# TP

IN=results/benchmarks/rtg_vcfeval/merged/tp_4callers.vcf.gz
OUT=results/benchmarks/rtg_vcfeval/merged/tp_4callers.csv

echo "CHROM,POS,REF,ALT,$(bcftools query -l "$IN" | paste -sd, -)" > "$OUT"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$IN" | tr '\t' ',' >> "$OUT"
```
