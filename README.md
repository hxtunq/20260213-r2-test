# Variant Calling Benchmarking Pipeline 

Pipeline so sánh hiệu suất 4 variant caller: GATK, DeepVariant, Strelka2, FreeBayes.

## Cập nhật chính của pipeline

- Pipeline mặc định chạy **toàn bộ chr22** (DeepVariant dùng `--model_type=WGS`; GATK/Strelka2/FreeBayes không còn giới hạn vùng exome ở bước gọi biến thể).
- Tự động phát hiện số CPU khả dụng và tự chọn `THREADS` (không cần chọn sẵn CPU trong script).
- Chỉ đo **runtime + CPU usage + Max RSS memory** ở 4 bước gọi biến thể (GATK/DeepVariant/Strelka2/FreeBayes), ghi vào `logs/resource_usage.tsv`.
- Giới hạn tài nguyên Docker theo cấu hình máy (`MAX_MEMORY=14G`) để tránh vượt khả năng máy.
- Thêm **Tầng B**: script đánh giá rủi ro chức năng của lỗi FP/FN và tính risk-weighted metrics.

## Quy trình thực hiện
...

Benchmark RTG vcfeval (Tầng A) sẽ được ghi tại:

- `results/benchmarks/rtg_vcfeval/<caller>/` (chi tiết TP/FP/FN)
- `results/benchmarks/rtg_vcfeval/vcfeval_summary.tsv` (bảng tổng hợp nhiều caller)

Kết quả Tầng B được ghi tại:

- `results/functional_risk/errors/<caller>/*_fp.vcf.gz`, `*_fn.vcf.gz`
- `results/functional_risk/risk_weighted_summary.tsv`
- `results/functional_risk/risk_weighted_details.tsv`

Tầng B sẽ so sánh theo 2 lớp callset:

- `raw` (PASS normalized)
- `wes_standard` (lọc kiểu WES thường dùng: `QUAL >= 30` và `DP >= 10` khi có trường DP)

và theo 5 mô hình risk-weight:

- `consequence`
- `alphamissense`
- `alphagenome`
- `varsage`
- `max_all` (lấy score lớn nhất trong các nguồn)

## Cấu trúc thư mục

```
.
├── config/
│   └── config.sh
├── scripts/
│   └── helper_functions.sh
├── data/
│   ├── reference/
│   │   ├── chr22.fa
│   │   └── chr22_non_N_regions.bed
│   └── simulated/
│       ├── *_truth.vcf.gz
│       ├── *_R1.fastq.gz
│       └── *_R2.fastq.gz
├── results/
│   ├── preprocessing/
│   └── variants/
│       ├── gatk/
│       ├── deepvariant/
│       ├── strelka2/
│       └── freebayes/
├── 02_preprocessing.sh
├── 03_variant_calling_gatk.sh
├── 04_variant_calling_deepvariant.sh
├── 05_variant_calling_strelka2.sh
├── 06_variant_calling_freebayes.sh
└── README.md
```

## Các tools được sử dụng:

- fastq, bwa, samtools, bcftools
- simuG, art_illumina
- gatk, freebayes, docker (cho DeepVariant và Strelka2)
