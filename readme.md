# Variant Benchmarking Repository

Scripts so sánh và đánh giá quy trình gọi biến thể của 4 công cụ gọi biến thể GATK, DeepVariant, Strelka2 và FreeBayes.

## 1. Cấu trúc thư mục ở giai đoạn ban đầu

```text
.
├── readme.md
├── workflow.md
├── stages/
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

## 2. Cấu trúc thư mục ở giai đoạn hoàn tất quy trình

```text
├── data/
│   ├── reference/          # reference genome, index, known-sites
│   └── simulated/          # truth VCF + FASTQ mô phỏng
├── log/                    # log dữ liệu chạy từ các tools
└── results/
    ├── preprocessing/      # BAM/metrics sau preprocessing
    |   └── fastqc_raw      # kiểm tra chất lượng 2 reads được giả lập từ ART
    ├── variants/           # VCF theo từng caller
    │   ├── gatk/
    │   ├── deepvariant/
    │   ├── strelka2/
    │   └── freebayes/
    ├── benchmarks/
    │   └── rtg_vcfeval/    # TP/FP/FN summary and merged
    └── functional_risk/    # risk-weighted summary/details
```

## 3. Thứ tự chạy pipeline

Xem tại file [workflow.md](https://github.com/hxtunq/20260213-r2-test/blob/main/workflow.md).




