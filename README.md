# Variant Calling Benchmarking Pipeline

Pipeline so sánh hiệu suất 4 variant caller: GATK, DeepVariant, Strelka2, FreeBayes

## Quy trình thực hiện

### PHẦN 1: Tạo cấu trúc folder

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

### PHẦN 2: Download và chuẩn bị dữ liệu

#### 2.1. Download Reference Genome (chr22 - hg38)

```bash
cd data/reference

# Download từ UCSC
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# Hoặc dùng curl
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz

# Giải nén
gunzip chr22.fa.gz

# Index reference
samtools faidx chr22.fa
bwa index chr22.fa
gatk CreateSequenceDictionary -R chr22.fa -O chr22.dict
```

#### 2.2. Tạo dữ liệu giả lập với Simutator

```bash
REF_FASTA="data/reference/chr22.fa"
SIM_DIR="data/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"

simutator mutate_fasta \
    --snps 7000 \
    --dels 2000:3 \
    --ins 2000:2 \
    --seed 42 \
    "${REF_FASTA}" \
    "${SIM_DIR}/${PREFIX}"
```

#### 2.3. Merge và chuẩn bị Truth VCF

```bash
REF_FASTA="data/reference/chr22.fa"
SIM_DIR="data/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
TRUTH_VCF="${SIM_DIR}/${PREFIX}_truth.vcf.gz"

# Merge tất cả VCF
bcftools concat ${SIM_DIR}/${PREFIX}*.original.vcf | bcftools sort -Oz -o ${TRUTH_VCF}
tabix -p vcf ${TRUTH_VCF}

# Tạo mutated FASTA
SAMPLE_NAME=$(bcftools query -l ${TRUTH_VCF} | head -n 1)
bcftools consensus -s "${SAMPLE_NAME}" -f "${REF_FASTA}" "${TRUTH_VCF}" > "${SIM_DIR}/${PREFIX}_mutated_combined.fa"
samtools faidx "${SIM_DIR}/${PREFIX}_mutated_combined.fa"

# Tạo file SNP và INDEL riêng
bcftools +fill-tags ${TRUTH_VCF} -Oz -o "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -- -t TYPE
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz"

bcftools view -i 'TYPE="snp"' "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -Oz -o "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
bcftools view -i 'TYPE="indel"' "${SIM_DIR}/${PREFIX}_truth_typed.vcf.gz" -Oz -o "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
```

#### 2.4. Tạo reads với ART Illumina

```bash
SIM_DIR="data/simulated"
PREFIX="SIMULATED_SAMPLE_chr22"
MUTATED_FASTA="${SIM_DIR}/${PREFIX}_mutated_combined.fa"

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

# Đổi tên và nén
mv "${SIM_DIR}/${PREFIX}_1.fq" "${SIM_DIR}/${PREFIX}_R1.fastq"
mv "${SIM_DIR}/${PREFIX}_2.fq" "${SIM_DIR}/${PREFIX}_R2.fastq"
gzip "${SIM_DIR}/${PREFIX}_R1.fastq"
gzip "${SIM_DIR}/${PREFIX}_R2.fastq"
```

#### 2.5. Tạo callable regions BED

```bash
REF_FAI="data/reference/chr22.fa.fai"
HIGH_CONF_BED="data/simulated/callable_regions.bed"

awk -v OFS='\t' '{print $1, 0, $2}' "${REF_FAI}" > "${HIGH_CONF_BED}"
```

### PHẦN 3: Chạy Pipeline

Sau khi hoàn thành PHẦN 1 và PHẦN 2, chạy các script theo thứ tự:

```bash
bash 02_preprocessing.sh
bash 03_variant_calling_gatk.sh
bash 04_variant_calling_deepvariant.sh
bash 05_variant_calling_strelka2.sh
bash 06_variant_calling_freebayes.sh
```

## Cấu trúc thư mục

```
.
├── config/
│   └── config.sh
├── scripts/
│   └── helper_functions.sh
├── data/
│   ├── reference/
│   │   └── chr22.fa
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

## Yêu cầu công cụ

- bwa
- samtools
- bcftools
- gatk
- fastp
- simutator
- art_illumina
- docker (cho DeepVariant và Strelka2)
- freebayes
