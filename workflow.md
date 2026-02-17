## Quy trình thực hiện

### Phần I. Tạo cấu trúc folder

```bash
git clone https://github.com/bitschif/variant-benchmarking.git
cd variant-benchmarking/
```

```bash
mkdir -p data/reference
mkdir -p data/simulated
mkdir -p results/preprocessing
mkdir -p results/variants/gatk
mkdir -p results/variants/deepvariant
mkdir -p results/variants/strelka2
mkdir -p results/variants/freebayes
mkdir -p logs
```

### Phần II. Download và chuẩn bị dữ liệu

#### 2.1. Download Reference Genome (chr22 - hg38)

```bash
#pwd: variant-benchmarking/data/reference

# download từ UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# giải nén
gunzip chr22.fa.gz

# index reference
samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

#### 2.2. Download Known Sites cho BQSR

Known sites giúp cải thiện chất lượng Base Quality Score Recalibration.

```bash
# pwd: variant-benchmarking/data/reference

# === dbSNP ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Extract chr22 và đảm bảo chromosome naming là "chr22"
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.vcf.gz

bcftools view -r chr22 -Oz -o dbsnp138.hg38.chr22.vcf.gz dbsnp138.hg38.vcf.gz
tabix -p vcf dbsnp138.hg38.chr22.vcf.gz

# === Mills and 1000G Gold Standard Indels ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -Oz -o Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz
tabix -p vcf Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz

# === 1000G Phase 1 SNPs ===
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# extract chr22
bcftools view -r chr22 1000G_phase1.snps.high_confidence.hg38.vcf.gz -Oz -o 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz
tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz

# xoá file không sử dụng đến
rm -f Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx
rm -f Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
rm -f 1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
```

**Một số lưu ý:**

```bash
# pwd = variant-benchmarking/data/reference

# reference genome từ UCSC sử dụng format "chr22"
# tất cả VCF files phải có tên chromosome khớp với tên ở file reference

# lệnh kiểm tra tên chromosome trong các file VCF
bcftools view -h file.vcf.gz | grep "^##contig"
tabix -l gile.vcf.gz # ở đồ án và repository này lấy tên là "chr22" làm chuẩn

# nếu gặp lỗi "chromosome not found", kiểm tra và đổi lại tên chromosome
# cách đổi tên "22" thành "chr22" nếu cần
bcftools annotate --rename-chrs chr_map.txt input.vcf.gz -Oz -o output.vcf.gz

# tạo file BED chứa vùng non-N của chromosome
python3 << 'EOF'
from Bio import SeqIO
import re

fasta = "chr22.fa"
output = "chr22_non_N_regions.bed"

with open(output, 'w') as out:
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq).upper()
        chrom = record.id
        for m in re.finditer(r'[ATGC]+', seq):
            start = m.start()
            end = m.end()
            if end - start >= 1000:
                out.write(f"{chrom}\t{start}\t{end}\n")

print("Done")
EOF
```

#### 2.3. Tạo đột biến SNPs và INDELs với simuG

```bash
# pwd: variant-benchmarking/data

# clone simuG về folder <data>
git clone https://github.com/yjx1217/simuG.git

# dùng simuG tạo đột biến
perl simuG/simuG.pl \
  -refseq reference/chr22.fa \
  -snp_count 60000 \
  -indel_count 30000 \
  -indel_min_len 1 \
  -indel_max_len 5 \
  -titv_ratio 2.0 \
  -prefix simulated/SIMULATED_SAMPLE_chr22
```

#### 2.4. Merge các file VCF output của simuG tạo truth VCF

```bash
# pwd: variant-benchmarking/data

# đảm bảo có .fai từ reference
samtools faidx reference/chr22.fa

# reheader từng file VCF với contig info
bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.vcf

bcftools reheader --fai reference/chr22.fa.fai \
  -o simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.vcf

# concat sau khi đã có header chuẩn
bcftools concat \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.SNP.reheader.vcf \
  simulated/SIMULATED_SAMPLE_chr22.refseq2simseq.INDEL.reheader.vcf | \
bcftools sort -Oz -o simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

tabix -p vcf simulated/SIMULATED_SAMPLE_chr22_truth.vcf.gz

```

#### 2.5. Tạo file fastq giả lập bằng ART Illumina

```bash
#pwd: variant-benchmarking/data

SIM_DIR="simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}.simseq.genome.fa"

art_illumina \
    -ss HS25 \
    -i "${MUTATED_FASTA}" \
    -p \
    -l 150 \
    -f 60 \
    -m 350 \
    -s 50 \
    -rs 42 \
    -o "${SIM_DIR}/${PREFIX}_" \
    -na

# đổi tên và nén
mv "${SIM_DIR}/${PREFIX}_1.fq" "${SIM_DIR}/${PREFIX}_R1.fastq"
mv "${SIM_DIR}/${PREFIX}_2.fq" "${SIM_DIR}/${PREFIX}_R2.fastq"
gzip "${SIM_DIR}/${PREFIX}_R1.fastq"
gzip "${SIM_DIR}/${PREFIX}_R2.fastq"
```

### Phần III. Tiền xử lý dữ liệu cho công đoạn gọi biến thể

Pipeline: FastQC → BWA-MEM → MarkDuplicates → BQSR

#### 3.1. FastQC - Raw reads quality control

```bash
mkdir -p results/preprocessing/fastqc_raw
fastqc -t 4 \
  -o results/preprocessing/fastqc_raw \
  data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  2>&1 | tee logs/fastqc_raw.log
```

#### 3.2. BWA-MEM - Alignment (directly from raw FASTQ after FastQC)

```bash
bwa mem \
  -t 4 \
  -R "@RG\tID:SIMULATED_SAMPLE\tSM:SIMULATED_SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  -M \
  data/reference/chr22.fa \
  data/simulated/SIMULATED_SAMPLE_chr22_R1.fastq.gz \
  data/simulated/SIMULATED_SAMPLE_chr22_R2.fastq.gz \
  2> logs/bwa_mem.log | \
  samtools sort -@ 4 -m 2G -o results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam
```

#### 3.3. GATK MarkDuplicates

```bash
gatk MarkDuplicates \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_aligned.bam \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  -M results/preprocessing/SIMULATED_SAMPLE_chr22_dup_metrics.txt \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX true \
  2>&1 | tee logs/markduplicates.log
```


#### 3.4. GATK BaseRecalibrator & ApplyBQSR

```bash
# Ensure known-sites VCFs match the "chr" naming convention if needed.
# If your VCFs use "22" instead of "chr22", run:
# scripts/rename_chromosomes.sh <input.vcf.gz> <output.vcf.gz>

gatk BaseRecalibrator \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --known-sites data/reference/dbsnp138.hg38.chr22.vcf.gz \
  --known-sites data/reference/Mills_and_1000G_gold_standard.indels.hg38.chr22.vcf.gz \
  --known-sites data/reference/1000G_phase1.snps.high_confidence.hg38.chr22.vcf.gz \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  2>&1 | tee logs/baserecalibrator.log

gatk ApplyBQSR \
  --java-options "-Xmx12G -XX:ParallelGCThreads=2" \
  -R data/reference/chr22.fa \
  -I results/preprocessing/SIMULATED_SAMPLE_chr22_marked.bam \
  --bqsr-recal-file results/preprocessing/SIMULATED_SAMPLE_chr22_recal.table \
  -O results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam \
  2>&1 | tee logs/applybqsr.log

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam
```

#### 3.5. Filter by mapping quality

```bash
samtools view \
  -@ 4 \
  -b \
  -q 20 \
  -F 1796 \
  results/preprocessing/SIMULATED_SAMPLE_chr22_recal.bam | \
  samtools sort -@ 4 -o results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam -

samtools index results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam
```

#### 3.6 Alignment statistics (samtools stats, mosdepth)

```bash
samtools stats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_stats.txt
samtools flagstat results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_flagstat.txt
samtools idxstats results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam \
  > results/preprocessing/SIMULATED_SAMPLE_chr22_idxstats.txt

# mosdepth for coverage
mosdepth -t 4 --by 1000 \
  results/preprocessing/SIMULATED_SAMPLE_chr22_coverage \
  results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam
```

#### 3.7 Export for downstream scripts

```bash
echo "FINAL_BAM=results/preprocessing/SIMULATED_SAMPLE_chr22_final.bam" > results/preprocessing/bam_path.sh
```

### Phần IV. Gọi biến thể

```bash
#pwd: variant-benchmarking

bash stages/03_variant_calling_gatk.sh
bash stages/04_variant_calling_deepvariant.sh
bash stages/05_variant_calling_strelka2.sh
bash stages/06_variant_calling_freebayes.sh
```

### Phần V. So sánh các công cụ gọi biến thể bằng RTG Tools VCFEval 

#### 5.1 Chạy VCFeval cho 4 file VCF output của 4 công cụ gọi biến thể sau khi gọi biến thể

```bash
# pwd: variant-benchmarking

# khai báo đường dẫn
BED=data/reference/chr22_non_N_regions.bed
REF=data/reference/chr22.fa
TRUTH_RAW=$(ls -1 data/simulated/*_truth.vcf.gz | head -n 1)

mkdir -p results/benchmarks/truth
rm -f results/benchmarks/truth/truth.gt.vcf.gz* results/benchmarks/truth/truth.gt.norm.vcf.gz*

# Tạo truth có FORMAT+GT hợp lệ (có khai báo ##FORMAT cho GT)
zcat "$TRUTH_RAW" | awk 'BEGIN{OFS="\t"; ff=0}
  /^##fileformat=/ {ff=1; print; next}
  /^##/ {print; next}
  /^#CHROM/ {
    if(ff==0) print "##fileformat=VCFv4.2"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print $0,"FORMAT","TRUTH"
    next
  }
  {print $0,"GT","1/1"}
' | bgzip -c > results/benchmarks/truth/truth.gt.vcf.gz

tabix -f -p vcf results/benchmarks/truth/truth.gt.vcf.gz

# Normalize lại truth (có GT) để dùng cho vcfeval
bcftools norm -f "$REF" -m -both results/benchmarks/truth/truth.gt.vcf.gz -Oz \
  -o results/benchmarks/truth/truth.gt.norm.vcf.gz

tabix -f -p vcf results/benchmarks/truth/truth.gt.norm.vcf.gz

# tạo SDF template cho RTG
rtg format -o results/benchmarks/ref/chr22.sdf "$REF"
# check
bcftools query -l results/benchmarks/truth/truth.gt.norm.vcf.gz
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/gatk/*_gatk_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/gatk \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/deepvariant/*_deepvariant_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/deepvariant \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/strelka2/*_strelka2_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/strelka2 \
  --threads 4
```

```bash
#pwd: variant-benchmarking
QUERY=$(ls -1 results/variants/freebayes/*_freebayes_pass.norm.vcf.gz | head -n 1)
RTG_MEM=14G rtg vcfeval \
  --baseline results/benchmarks/truth/truth.gt.norm.vcf.gz \
  --calls "$QUERY" \
  --template results/benchmarks/ref/chr22.sdf \
  --bed-regions "$BED" \
  --output results/benchmarks/rtg_vcfeval/freebayes \
  --threads 4
```

#### 5.2 Hợp nhất các file output fn/fp/tp sau khi dùng rtgtool/vcfeval của 4 công cụ gọi biến thể

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

### 5.3 Chuyển file VCF sang CSV để phục vụ cho công đoạn xử lý dữ liệu và trực quan hoá

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
