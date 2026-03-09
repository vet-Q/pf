# Plasmodium falciparum (Pf) 분석 파이프라인

Nanopore 시퀀싱으로 생성된 _Plasmodium falciparum_ FASTQ 데이터를 정렬(alignment)부터 변이 호출(variant calling), 합의 서열(consensus), 커버리지 분석, 시각화까지 한 번에 처리하는 **Nextflow 기반 파이프라인**입니다.

---

## 📁 프로젝트 구조

```
pf/
├── main.nf                  ← 메인 파이프라인 스크립트
├── nextflow.config          ← 모든 파라미터 기본값 + 실행 프로파일
├── run.sh                   ← 한 줄 실행 스크립트
├── setup.sh                 ← 환경 점검 스크립트
│
├── configs/
│   └── params.yaml          ← (참고용) 파라미터 설명 파일
│
├── data/                    ← 📂 입력 데이터 (barcode별 FASTQ)
│   ├── barcode01/
│   │   ├── *.fastq.gz
│   │   └── ...
│   ├── barcode02/
│   └── ... (barcode24까지)
│
├── resources/               ← 참조 파일
│   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta   ← 레퍼런스 게놈
│   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta.fai
│   └── amplicons.bed        ← 유전자 좌표 (CRT, DHFR, ...)
│
├── modules/                 ← Nextflow 모듈 (각 분석 단계)
│   ├── minimap2_align.nf    ← FASTQ → BAM 정렬
│   ├── pf_depth.nf          ← 깊이 분석
│   ├── ivar_consensus.nf    ← 합의 서열 생성
│   ├── call_variants.nf     ← bcftools 변이 호출
│   ├── gene_bins_nf.nf      ← 유전자별 read-length bin 분석
│   ├── gene_bins_plot_nf.nf ← 멀티패널 플롯
│   ├── gene_1to8_plot.nf    ← 1×8 유전자 플롯
│   ├── af_dp_figures.nf     ← AF/DP 산점도 피규어
│   └── fastq_qc.nf          ← FastQC 품질 관리
│
├── bin/                     ← 헬퍼 스크립트 (Shell + R)
│   ├── call_bcftools.sh
│   ├── call_ivar.sh
│   ├── make_consensus_bcftools.sh
│   ├── make_consensus_ivar.sh
│   ├── make_af_dp_figures.R
│   ├── plot_depth_per_gene.R
│   ├── plot_gene_multipanel.R
│   └── plot_gene_1x8.R
│
└── results/                 ← 📂 출력 결과 (자동 생성)
    ├── alignment/minimap2/  ← 정렬된 BAM 파일
    ├── pf_depth/            ← 깊이 분석 결과
    ├── calling/ivar/        ← iVar 합의 서열 + 변이
    ├── calling/bcftools/    ← bcftools VCF 파일
    ├── gene_bins/           ← 유전자별 bin 분석
    ├── figures/afdp/        ← AF/DP 피규어
    └── pipeline_info/       ← Nextflow 실행 리포트
```

---

## 🔧 사전 요구사항

### 필수 소프트웨어

| 도구 | 용도 | 설치 확인 |
|------|------|-----------|
| **Java** 11+ | Nextflow 실행 | `java -version` |
| **Nextflow** ≥ 22.0 | 파이프라인 엔진 | `nextflow -version` |
| **samtools** | BAM 처리 / 인덱싱 | `samtools --version` |
| **minimap2** | FASTQ → BAM 정렬 | `minimap2 --version` |
| **bcftools** | 변이 호출 (VCF) | `bcftools --version` |
| **ivar** | 합의 서열 생성 | `ivar version` |
| **Rscript** ≥ 4.0 | 시각화 (ggplot2 등) | `Rscript --version` |

### 설치 방법 (Ubuntu/Debian)

```bash
# Java
sudo apt install default-jdk

# Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# 생물정보학 도구 (conda 추천)
conda install -c bioconda samtools minimap2 bcftools ivar

# R 패키지 (R 콘솔에서)
# install.packages(c("ggplot2", "dplyr", "patchwork", "stringr", "purrr", "optparse"))
# BiocManager::install("vcfR")
```

### 환경 점검

```bash
./setup.sh
```

이 명령으로 모든 도구가 설치되어 있는지 한눈에 확인할 수 있습니다.

---

## 🚀 실행 방법

### 방법 1: run.sh 스크립트 (가장 간단)

```bash
cd /mnt/disk1/NGS/Plasmodium_falciparum/analysis/pf
./run.sh
```

### 방법 2: Nextflow 직접 실행

```bash
cd /mnt/disk1/NGS/Plasmodium_falciparum/analysis/pf
nextflow run main.nf -profile standard
```

### 방법 3: 파라미터 오버라이드

```bash
# 예) 정렬만 실행
nextflow run main.nf -profile standard \
  --do_alignment true \
  --do_depth false \
  --do_ivar false \
  --do_variants false \
  --do_afdp false

# 예) 이미 정렬된 BAM에서 시작
nextflow run main.nf -profile standard \
  --do_alignment false \
  --bam_root /path/to/bam/files \
  --bam_glob "*.sorted.bam"

# 예) 특정 유전자만 분석
nextflow run main.nf -profile standard \
  --genes "CRT K13 MDR1"
```

---

## ⚙️ 파이프라인 단계 설명

파이프라인은 아래 순서로 실행됩니다. 각 단계는 `nextflow.config`에서 on/off 가능합니다.

```
FASTQ 입력 (data/barcode*/*.fastq.gz)
    │
    ▼
┌─────────────────────────────────────┐
│ ① MINIMAP2_ALIGN (do_alignment)     │  FASTQ → 정렬된 BAM
│    preset: map-ont (Nanopore)       │
└─────────────────────────────────────┘
    │
    ▼  BAM 파일
    ├──────────────────────┬───────────────────┬─────────────────────┐
    ▼                      ▼                   ▼                     ▼
┌──────────┐       ┌──────────────┐   ┌──────────────┐     ┌──────────────┐
│ ② 깊이   │       │ ③ iVar       │   │ ④ 유전자별   │     │ ⑤ bcftools   │
│  분석     │       │  합의서열    │   │  bin 분석    │     │  변이호출    │
│(do_depth)│       │ (do_ivar)    │   │(do_gene_bins)│     │(do_variants) │
└──────────┘       └──────────────┘   └──────┬───────┘     └──────┬───────┘
                                              │                    │
                                              ▼                    ▼
                                     ┌──────────────┐     ┌──────────────┐
                                     │ 멀티패널 플롯 │     │ ⑥ AF/DP     │
                                     │(do_multipanel)│     │  피규어     │
                                     └──────────────┘     │ (do_afdp)   │
                                                          └──────────────┘
```

### 각 단계 상세

| 단계 | 파라미터 | 설명 |
|------|----------|------|
| **① 정렬** | `do_alignment: true` | minimap2로 FASTQ를 레퍼런스에 정렬하여 BAM 생성 |
| **② 깊이 분석** | `do_depth: true` | amplicon 영역의 position별 시퀀싱 깊이를 계산하고 유전자별 플롯 생성 |
| **③ iVar 합의서열** | `do_ivar: true` | 유전자별 합의 서열(consensus FASTA)과 변이(TSV) 생성 |
| **④ 유전자 bin 분석** | `do_gene_bins: true` | Read 길이를 0-1k, 1-2k, 2-3k, 3-4k, 4k+ 구간으로 분류하여 깊이 분석 |
| **⑤ 변이 호출** | `do_variants: true` | bcftools로 SNP 호출 → VCF.gz 생성 + bcftools consensus |
| **⑥ AF/DP 피규어** | `do_afdp: true` | 모든 유전자의 Allele Frequency vs Depth 산점도/히스토그램 생성 |

---

## ⚙️ 주요 파라미터 (nextflow.config)

모든 파라미터는 `nextflow.config`의 `params {}` 블록에 정의되어 있습니다.

### 경로 설정

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `fastq_root` | `${projectDir}/data` | FASTQ 파일이 있는 루트 디렉토리 |
| `fastq_glob` | `barcode*/*.fastq.gz` | FASTQ 파일 패턴 |
| `bam_root` | `${projectDir}/results/alignment/minimap2` | BAM 파일 경로 (정렬 건너뛸 때) |
| `outdir` | `${projectDir}/results` | 결과 출력 디렉토리 |

### 분석 도구 설정

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `minimap2_preset` | `map-ont` | 정렬 프리셋 (`map-ont`=Nanopore, `sr`=short reads, `map-pb`=PacBio) |
| `genes` | `CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII` | 분석 대상 약제내성 유전자 |
| `ivar_threshold` | `0.6` | iVar 합의서열 생성 빈도 임계값 |
| `snp_only` | `true` | SNP만 호출할지 여부 |
| `max_depth` | `10000` | bcftools mpileup 최대 깊이 |

### 파이프라인 토글 (on/off)

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `do_alignment` | `true` | FASTQ → BAM 정렬 수행 |
| `do_depth` | `true` | 깊이 분석 |
| `do_ivar` | `true` | iVar 합의서열 생성 |
| `do_gene_bins` | `true` | 유전자별 read-length bin 분석 |
| `do_variants` | `true` | bcftools 변이 호출 |
| `do_afdp` | `true` | AF/DP 피규어 생성 |
| `do_fastqc` | `false` | FastQC 품질 보고서 |
| `do_multipanel` | `false` | 멀티패널 플롯 |
| `do_gene_1x8` | `false` | 1×8 유전자 플롯 |

---

## 📊 출력 결과

```
results/
├── alignment/minimap2/
│   ├── barcode01.Pf3D7.sorted.bam     ← 정렬된 BAM
│   ├── barcode02.Pf3D7.sorted.bam
│   └── ...
│
├── pf_depth/
│   └── barcode01/
│       ├── all_amplicons.depth.txt     ← 전체 amplicon 깊이
│       ├── depth_by_gene/              ← 유전자별 깊이 텍스트
│       │   ├── CRT_depth.txt
│       │   └── ...
│       └── plots_by_gene/              ← 유전자별 깊이 플롯 (PDF)
│           ├── CRT_depth.pdf
│           └── ...
│
├── calling/ivar/
│   └── barcode01/
│       └── CRT/
│           ├── barcode01.CRT.ivar.consensus.fasta  ← 합의 서열
│           ├── barcode01_CRT_ivar.all.tsv           ← 모든 변이
│           ├── barcode01_CRT_ivar.snp.tsv           ← SNP만
│           └── depth.tsv                            ← position별 깊이
│
├── calling/bcftools/
│   └── barcode01/
│       └── CRT/
│           ├── barcode01.CRT.bcftools.vcf.gz       ← VCF 파일
│           ├── barcode01.CRT.bcftools.consensus.fasta
│           └── region.bed
│
├── gene_bins/
│   └── barcode01/
│       └── CRT/
│           ├── bin_0-1k.bam           ← 길이별 BAM
│           ├── bin_0-1k_depth.txt     ← 길이별 깊이
│           └── readlen.txt            ← read 길이 분포
│
├── figures/afdp/
│   ├── Figure1_by_gene.pdf            ← 유전자별 AF vs DP 산점도
│   └── Figure2_overall.pdf            ← 전체 AF vs DP 산점도
│
└── pipeline_info/
    ├── report.html                    ← Nextflow 실행 리포트
    └── timeline.html                  ← 실행 시간 타임라인
```

---

## 🔄 재실행 / 이어서 실행

```bash
# 캐시된 결과 활용하여 이어서 실행 (실패한 단계부터)
nextflow run main.nf -profile standard -resume

# 완전히 새로 실행 (모든 캐시 삭제)
rm -rf work/ .nextflow/ .nextflow.log*
nextflow run main.nf -profile standard
```

---

## 🐛 문제 해결

### "No FASTQ files found" 오류
```bash
# data 폴더에 파일이 있는지 확인
ls data/barcode01/*.fastq.gz | head -3
```

### "Command not found" 오류
```bash
# 필요한 도구가 설치되어 있는지 확인
./setup.sh
```

### 메모리 부족
```bash
# nextflow.config에서 리소스 조정
# params.max_memory = '64.GB'
# 또는 실행 시 오버라이드
nextflow run main.nf -profile standard --max_memory '64.GB'
```

### 특정 단계만 다시 실행
```bash
# 예) 변이 호출만 다시 실행 (정렬은 캐시 사용)
nextflow run main.nf -profile standard -resume --do_variants true
```

---

## 📜 라이선스

연구 목적으로 자유롭게 사용 가능합니다.
