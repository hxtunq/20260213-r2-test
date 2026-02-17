# Variant Benchmarking Repository

Repository này đã chạy hoàn chỉnh và được tổ chức lại để:
- Giữ các tài liệu/manual quan trọng ngay ở thư mục gốc.
- Gom toàn bộ script chạy pipeline vào một nơi duy nhất để dễ thao tác.
- Trình bày dữ liệu/kết quả rõ ràng để người ngoài xem nhanh.

## 1) Cấu trúc thư mục sau khi chuẩn hoá

```text
.
├── 02_preprocessing_manual.sh
├── 08_vcf_merge_manual.md
├── README.md
├── workflow.md
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

## 2) Quy ước đặt chỗ dữ liệu để public dễ xem

Vì bạn public toàn bộ data, nên tách rõ phần "input" và "output" như sau:

```text
data/
├── reference/          # reference genome, index, known-sites
├── simulated/          # truth VCF + FASTQ mô phỏng
└── metadata/           # mô tả sample, checksum, phiên bản tool

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

## 3) Gợi ý phân nhánh để người ngoài xem kết quả nhanh

- `main`: mã nguồn + script + tài liệu (không đẩy file output nặng).
- `data-public`: toàn bộ dữ liệu đầu vào đã public (`data/`).
- `results-public`: toàn bộ output đã chạy (`results/`, `logs/`).

Nếu dữ liệu lớn, nên dùng Git LFS cho `*.bam`, `*.vcf.gz`, `*.tbi`, `*.sdf`.

## 4) Thứ tự chạy pipeline hiện tại

Từ thư mục gốc repo:

```bash
bash 02_preprocessing_manual.sh
bash pipeline/03_variant_calling_gatk.sh
bash pipeline/04_variant_calling_deepvariant.sh
bash pipeline/05_variant_calling_strelka2.sh
bash pipeline/06_variant_calling_freebayes.sh
bash pipeline/07_benchmark_rtg_vcfeval.sh
```

## 5) File giữ ở thư mục gốc (theo yêu cầu)

- `02_preprocessing_manual.sh`
- `08_vcf_merge_manual.md`
- `README.md`
- `workflow.md`

