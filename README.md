# Variant Benchmarking Repository

Scripts so sánh và đánh giá quy trình gọi biến thể của 4 công cụ gọi biến thể GATK, DeepVariant, Strelka2 và FreeBayes.

## 1. Cấu trúc thư mục tại trạng thái ban đầu

```text
.
├── readme.md
├── workflow.md
├── data/
│   ├── reference/          # reference genome, index, known-sites
│   └── simulated/          # truth VCF + FASTQ mô phỏng
├── pipeline/
│   ├── 03_variant_calling_gatk.sh
│   ├── 04_variant_calling_deepvariant.sh
│   ├── 05_variant_calling_strelka2.sh
│   ├── 06_variant_calling_freebayes.sh
│   └── 07_benchmark_rtg_vcfeval.sh
├── config/
│   └── config.sh
└── scripts/
    ├── helper_functions.sh
    ├── normalize_vcf.sh
    ├── rename_chromosomes.sh
    └── risk_weighted_eval.py
```

## 2. Cấu trúc thư mục sau khi chạy lệnh

```text
results/
├── preprocessing/      # BAM/metrics sau preprocessing
├── variants/           # VCF theo từng caller
│   ├── gatk/
│   ├── deepvariant/
│   ├── strelka2/
│   └── freebayes/
├── benchmarks/
│   └── rtg_vcfeval/    # TP/FP/FN + summary
└── functional_risk/    # risk-weighted summary/details
```

## 3. Thứ tự chạy pipeline hiện tại

Từ thư mục gốc repo:

```bash
bash pipeline/03_variant_calling_gatk.sh
bash pipeline/04_variant_calling_deepvariant.sh
bash pipeline/05_variant_calling_strelka2.sh
bash pipeline/06_variant_calling_freebayes.sh
bash pipeline/07_benchmark_rtg_vcfeval.sh
```

