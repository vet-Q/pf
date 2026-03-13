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
│                              ※ 실행 시 -params-file 옵션으로 지정해야 적용됨
│
├── data/                    ← 📂 입력 데이터 (barcode별 FASTQ)
│   ├── barcode01/
│   │   └── *.fastq.gz
│   └── ... (barcode24까지)
│
├── resources/               ← 참조 파일
│   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta   ← 레퍼런스 게놈
│   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta.fai
│   ├── amplicons.bed        ← amplicon 영역 좌표 (8개 유전자)
│   └── cds_coords.tsv       ← CDS exon 좌표 (번역용, PlasmoDB-55 기반)
│
├── modules/                 ← Nextflow 모듈 (각 분석 단계)
│   ├── minimap2_align.nf      ← FASTQ → BAM 정렬
│   ├── amplicon_depth.nf      ← 깊이 분석
│   ├── ivar_consensus.nf      ← 합의 서열 생성
│   ├── call_variants.nf       ← bcftools 변이 호출
│   ├── collect_variants.nf    ← 전체 샘플 VCF 병합 → variant_AF_all.tsv
│   ├── splice_translate.nf    ← CDS 접합 + 번역 (샘플별)
│   ├── align_variants.nf      ← MAFFT 정렬 + 변이표 생성 (유전자별)
│   ├── gene_bin_depth.nf      ← 유전자별 read-length bin 분석
│   ├── gene_multipanel_plot.nf← 멀티패널 플롯
│   ├── gene_1x8_plot.nf       ← 1×8 유전자 플롯
│   ├── af_dp_figures.nf       ← AF/DP 산점도 피규어
│   ├── resistance_table.nf    ← 약제내성 마커 요약표
│   └── fastq_qc.nf            ← FastQC 품질 관리
│
├── bin/                     ← 헬퍼 스크립트 (Shell + Python + R)
│   ├── call_bcftools.sh
│   ├── call_ivar.sh
│   ├── make_consensus_bcftools.sh
│   ├── make_consensus_ivar.sh
│   ├── splice_and_translate.py    ← CDS exon 추출 + 접합 + 번역 (Python)
│   ├── gene_msa_variants.R        ← 유전자별 MSA + 변이표 생성 (R)
│   ├── collect_variants.R         ← 전체 샘플 VCF 병합 → variant_AF_all.tsv (R)
│   ├── make_resistance_table.R    ← 약제내성 마커 요약표 생성 (R)
│   ├── make_af_dp_figures.R
│   ├── plot_depth_per_gene.R
│   ├── plot_gene_multipanel.R
│   └── plot_gene_1x8.R
│
└── results/                 ← 📂 출력 결과 (자동 생성)
    ├── alignment/minimap2/  ← 정렬된 BAM 파일
    ├── pf_depth/            ← 깊이 분석 결과
    ├── calling/ivar/        ← iVar 합의 서열 + 변이
    ├── calling/bcftools/    ← bcftools VCF + consensus FASTA
    ├── calling/variant_AF_all.tsv ← 전체 샘플 병합 변이표
    ├── translation/         ← CDS 번역 + MAFFT 정렬 + 변이표
    ├── gene_bins/           ← 유전자별 bin 분석
    ├── figures/afdp/        ← AF/DP 피규어
    ├── resistance/          ← 약제내성 마커 요약표
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
| **mafft** | 다중 서열 정렬 (번역 기능) | `mafft --version` |
| **python3** ≥ 3.8 | CDS 접합·번역·변이표 스크립트 | `python3 --version` |
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
conda install -c bioconda samtools minimap2 bcftools ivar mafft

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
                                              ▼                    │
                                     ┌──────────────┐             │
                                     │ 멀티패널 플롯 │             │
                                     │(do_multipanel)│             │
                                     └──────────────┘             │
                                                                   ├──────────────────┐
                                                                   ▼                  ▼
                                                          ┌──────────────┐  ┌─────────────────────┐
                                                          │ ⑥ AF/DP     │  │ ⑦ CDS 번역/변이표   │
                                                          │  피규어     │  │  (do_translate)     │
                                                          │ (do_afdp)   │  │                     │
                                                          └──────────────┘  │ SPLICE_TRANSLATE    │
                                                                            │  ↓ (샘플별)         │
                                                                            │ ALIGN_VARIANTS      │
                                                                            │  ↓ (유전자별)       │
                                                                            │ NT/AA MSA + 변이표  │
                                                                            └─────────────────────┘
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
| **⑦ CDS 번역·변이표** | `do_translate: false` | consensus에서 exon만 추출→접합→번역→MAFFT 정렬→변이표 생성 (아래 상세 참조) |
| **⑧ 약제내성 마커표** | `do_resistance: false` | 전체 VCF를 병합 후 알려진 내성 SNP 위치만 추출하여 샘플×마커 요약표 생성 |

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
| `do_translate` | `false` | CDS 번역·다중 정렬·변이표 생성 (⑦단계, mafft 필요) |
| `do_resistance` | `false` | 알려진 약제내성 마커 요약표 생성 (⑧단계) |
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
├── translation/
│   ├── per_sample/                    ← 샘플별 spliced 서열
│   │   ├── barcode01.CRT.nt.fasta    ← exon 접합 CDS (핵산)
│   │   ├── barcode01.CRT.aa.fasta    ← 번역 단백질 서열
│   │   ├── barcode01.CRT.af.tsv      ← IUPAC 해소 AF 정보표
│   │   ├── barcode01.DHFR.nt.fasta
│   │   ├── barcode01.DHFR.aa.fasta
│   │   ├── barcode01.DHFR.af.tsv
│   │   └── ...
│   ├── CRT/
│   │   ├── CRT.nt.combined.fasta     ← 3D7 + 전체 샘플 NT 합본
│   │   ├── CRT.nt.msa.fasta          ← MAFFT 정렬 NT
│   │   ├── CRT.nt.variants.tsv       ← 핵산 변이표
│   │   ├── CRT.aa.combined.fasta     ← 3D7 + 전체 샘플 AA 합본
│   │   ├── CRT.aa.msa.fasta          ← MAFFT 정렬 AA
│   │   └── CRT.aa.variants.tsv       ← 아미노산 변이표 (★ 핵심 결과)
│   ├── DHFR/
│   │   └── ...
│   └── ...
│
├── resistance/
│   └── resistance_markers.tsv         ← 샘플 × 마커 약제내성 요약표
│
├── calling/
│   └── variant_AF_all.tsv             ← 전체 샘플 병합 변이표
│
└── pipeline_info/
    ├── report.html                    ← Nextflow 실행 리포트
    └── timeline.html                  ← 실행 시간 타임라인
```

---

## 💊 약제내성 마커 요약표 (`do_resistance`)

`do_resistance: true`로 활성화하면 bcftools VCF 결과에서 알려진 약제내성 SNP 위치만 추출하여
**샘플 × 마커**의 wide-format TSV를 자동 생성합니다.

### 처리 흐름

```
CALL_VARIANTS → 샘플별/유전자별 .vcf.gz
       ↓  (collect_variants.R)
COLLECT_VARIANTS → calling/variant_AF_all.tsv  (전체 샘플 병합)
       ↓  (make_resistance_table.R)
RESISTANCE_TABLE → resistance/resistance_markers.tsv
```

### 출력 예시 (`resistance_markers.tsv`)

| barcode | CRT_C72S | CRT_V73I | CRT_M74I | CRT_N75E | CRT_K76T | DHFR_S108N | K13_C580Y | … |
|---------|----------|----------|----------|----------|----------|------------|-----------|---|
| barcode01 | C | V | M | N | T(1.00) | N(0.99) | C | … |
| barcode02 | C | V | I(0.95) | E(0.95) | T(1.00) | N(0.98) | Y(0.93) | … |
| barcode03 | C | V | M | N | K | S(0.82) | C | … |

- 레퍼런스(3D7)와 동일 → ref 아미노산 단일 문자 (예: `K`, `C`)
- 명확한 변이 (AF ≥ 0.5) → `alt_aa(AF)` (예: `T(1.00)`, `Y(0.93)`)
- 혼합 (0.2 ≤ AF < 0.5) → `ref/alt(AF)` (예: `K/T(0.35)`)
- 노이즈 (AF < 0.2) → ref 아미노산

### 수록 마커 목록

#### CRT (PF3D7_0709000) — Chloroquine / Piperaquine

참조 하플로타입 (3D7): **C**72-**V**73-**M**74-**N**75-**K**76 (CVMNK) → 내성: CVIET

| 마커 | aa 위치 | ref → alt | 약제 | 비고 |
|------|---------|-----------|------|---------|
| CRT_C72S | 72 | C → S | CQ | 일부 동남아 분리주에서 보고 |
| CRT_V73I | 73 | V → I | CQ | 드물게 보고 |
| CRT_M74I | 74 | M → I | CQ | CVIET 하플로타입 |
| CRT_N75E | 75 | N → E | CQ | CVIET 하플로타입 |
| **CRT_K76T** | **76** | **K → T** | **CQ** | **핵심 CQ 내성 마커** |
| CRT_A220S | 220 | A → S | CQ | 아프리카 내성 분리주 |
| CRT_Q271E | 271 | Q → E | CQ | |
| CRT_N326S | 326 | N → S | CQ | |
| CRT_I356T | 356 | I → T | CQ | |
| CRT_R371I | 371 | R → I | CQ | |

#### DHFR (PF3D7_0417200) — Pyrimethamine

| 마커 | aa 위치 | ref → alt | 약제 | 비고 |
|------|---------|-----------|------|---------|
| DHFR_N51I | 51 | N → I | PYR | 이중 변이 (51+108) |
| DHFR_C59R | 59 | C → R | PYR | 삼중 변이 |
| **DHFR_S108N** | **108** | **S → N** | **PYR** | **핵심 PYR 내성 마커** |
| DHFR_I164L | 164 | I → L | PYR | 사중 변이 (고도내성) |

#### DHPS (PF3D7_0810800) — Sulfadoxine

| 마커 | aa 위치 | ref → alt | 약제 | 비고 |
|------|---------|-----------|------|---------|
| DHPS_S436A | 436 | S → A | SDX | |
| **DHPS_A437G** | **437** | **A → G** | **SDX** | **핵심 SDX 내성 마커** |
| DHPS_K540E | 540 | K → E | SDX | 삼중 변이 (434+437+540) |
| DHPS_A581G | 581 | A → G | SDX | 오중 변이 (고도내성) |
| DHPS_A613S | 613 | A → S | SDX | 오중 변이 |

#### K13 (PF3D7_1343700) — Artemisinin

WHO validated markers + candidate markers 포함

| 마커 | aa 위치 | ref → alt | 약제 | 비고 |
|------|---------|-----------|------|---------|
| K13_F446I | 446 | F → I | ART | WHO validated |
| K13_M476I | 476 | M → I | ART | WHO validated |
| K13_Y493H | 493 | Y → H | ART | WHO validated |
| K13_G533V | 533 | G → V | ART | WHO validated |
| K13_R539T | 539 | R → T | ART | WHO validated |
| K13_I543T | 543 | I → T | ART | WHO validated |
| K13_P553L | 553 | P → L | ART | Candidate |
| K13_R561H | 561 | R → H | ART | WHO validated (동남아 주요) |
| K13_P574L | 574 | P → L | ART | Candidate |
| **K13_C580Y** | **580** | **C → Y** | **ART** | **WHO validated (동아프리카 주요)** |
| K13_A675V | 675 | A → V | ART | WHO validated |

#### MDR1 (PF3D7_0523000) — Lumefantrine / Amodiaquine

| 마커 | aa 위치 | ref → alt | 약제 | 비고 |
|------|---------|-----------|------|---------|
| **MDR1_N86Y** | **86** | **N → Y** | **LUM/AQ** | **lumefantrine 내성 관련** |
| MDR1_Y184F | 184 | Y → F | LUM/AQ | |
| MDR1_S1034C | 1034 | S → C | LUM/AQ | |
| MDR1_N1042D | 1042 | N → D | LUM/AQ | |
| MDR1_D1246Y | 1246 | D → Y | LUM/AQ | |

### 실행 방법

```bash
# Variant calling과 함께 한 번에 실행
nextflow run main.nf -profile standard \
  --do_alignment false \
  --do_variants true \
  --do_resistance true

# 기존 VCF 결과에서 resistance table만 추가 실행
nextflow run main.nf -profile standard \
  --do_alignment false \
  --do_depth false \
  --do_ivar false \
  --do_gene_bins false \
  --do_afdp false \
  --do_variants false \
  --do_resistance true
```

> `do_variants=false`이더라도 `results/calling/bcftools/` 폴더가 존재하면
> 자동으로 기존 VCF를 읽어 `COLLECT_VARIANTS → RESISTANCE_TABLE`을 실행합니다.

---

## 🧬 CDS 번역 및 변이표 생성 (`do_translate`)

### 배경: amplicon BED ≠ ORF 범위

`amplicons.bed`에 정의된 amplicon 범위는 amplicon 프라이머 설계 기준이므로
실제 ORF(open reading frame)보다 넓습니다. 또한 일부 유전자(특히 **CRT**)는
genomic 서열 안에 여러 개의 intron이 존재하여, BED 범위의 consensus를 그대로
번역하면 frame이 깨집니다.

이 단계는 다음 문제를 모두 자동으로 처리합니다:

| 문제 | 처리 방식 |
|---|---|
| BED 범위 > ORF | `cds_coords.tsv`의 exon 좌표로 ORF만 잘라냄 |
| intron 포함 | genomic 좌표로 exon 구간만 추출하여 접합 |
| minus strand 유전자 | exon을 역방향으로 접합 후 역상보(RC) |
| 다중 exon (예: CRT 13개) | exon을 transcript 순서(5'→3')로 이어 붙임 |

### 대상 유전자 exon 구조

| 유전자 | PlasmoDB ID | exon 수 | strand | 단백질 길이 |
|---|---|---|---|---|
| **CRT** | PF3D7_0709000 | 13 | + | 424 aa |
| **DHFR** | PF3D7_0417200 | 1 | + | 608 aa |
| **DHPS** | PF3D7_0810800 | 3 | + | 706 aa |
| **K13** | PF3D7_1343700 | 1 | − | 726 aa |
| **MDR1** | PF3D7_0523000 | 1 | + | 1419 aa |
| **MSP2** | PF3D7_0206800 | 1 | − | 272 aa |
| **PMI** | PF3D7_1407900 | 1 | + | 452 aa |
| **PMIII** | PF3D7_1408100 | 1 | + | 451 aa |

> CDS exon 좌표는 `resources/cds_coords.tsv`에 0-based BED 형식으로 저장되어 있습니다.
> PlasmoDB-55 GFF annotation 기준으로 추출하였으며 PlasmoDB-67 게놈과 호환됩니다.

### 처리 흐름

```
[bcftools consensus FASTA]
  (BED 전체 범위, intron 포함)
          │
          │  splice_and_translate.py  (Python)
          │  1) FASTA 헤더에서 genomic offset 계산
          │  2) cds_coords.tsv에서 exon 좌표 읽기
          │  3) exon 구간만 잘라냄 (intron 제거)
          │  4) + strand: 좌→우 접합
          │     - strand: 우→좌 접합 후 역상보(RC)
          │  5) 표준 유전부호로 번역 (N→X 처리)
          ▼
  barcode01.CRT.nt.fasta   ← spliced CDS (핵산)
  barcode01.CRT.aa.fasta   ← 번역 단백질

[3D7 reference 추출]
  splice_and_translate.py --ref-fasta  (Python)
  → samtools faidx로 각 exon 구간 추출 후 접합·번역
  → 3D7 기준 서열 생성

[전체 샘플 + 3D7 합본 → MAFFT 정렬 → 변이표]
  gene_msa_variants.R  (R)
  → 유전자별로 모든 샘플 FASTA 수집
  → 3D7 reference 을 첫 번째로 붙여 combined FASTA 생성
  → mafft --auto 로 NT / AA 각각 정렬
  → 3D7 대비 variant 위치를 TSV 로 저장

최종 출력:
  CRT.aa.variants.tsv  ← pos / ref / barcode01 / barcode02 / ...
  barcode01.CRT.af.tsv ← CDS 내 IUPAC 위치별 AF 정보 (아래 참조)
```

### 변이표 예시 (CRT.aa.variants.tsv)

```
gene  pos  ref  barcode01  barcode02  barcode03  ...
CRT   72   C    C          S          S
CRT   74   M    M          M          I
CRT   75   N    N          E          E
CRT   76   K    T          T          T          ← K76T (CQ 내성 마커)
CRT   220  A    A          S          S          ← A220S
```

### 실행 방법

```bash
# nextflow.config에서 활성화 (do_variants도 true여야 함)
nextflow run main.nf -profile standard \
  --do_variants true \
  --do_translate true

# 또는 커맨드라인 직접 지정
nextflow run main.nf -profile standard --do_translate true

# do_variants 없이 기존 calling 결과에서 단독 실행 (-resume으로 CALL_VARIANTS 캐시 재사용)
nextflow run main.nf -profile standard \
  --do_translate true \
  --do_alignment false \
  --do_depth false \
  --do_ivar false \
  --do_gene_bins false \
  --do_variants false \
  --do_afdp false \
  -resume
```

> `do_variants=false`이더라도 `results/calling/bcftools/` 폴더가 존재하면
> 자동으로 기존 결과를 읽어 SPLICE_TRANSLATE → ALIGN_VARIANTS를 실행합니다.

---

### 🔬 IUPAC 염기 해소 (AF 기반)

#### 배경

`bcftools consensus`는 heterozygous site(REF와 ALT 둘 다 존재)를 **IUPAC 모호성 코드**로 표기합니다.
예: Y = C/T 혼재, R = A/G 혼재, K = G/T 혼재

이 경우 코돈 안에 IUPAC 코드가 있으면 번역 시 **X**(unknown)가 되어 아미노산 변이 분석이 불가능합니다.
또한 단순한 heterozygous call인지, **혼합감염(mixed infection)**인지 구분이 필요합니다.

#### 해소 방법

`splice_and_translate.py`는 bcftools VCF의 `FORMAT/AD`(allele depth)를 읽어
**실제 allele frequency(AF)**를 계산한 뒤 아래 규칙으로 IUPAC를 해소합니다:

| 조건 | 결정 | 의미 |
|---|---|---|
| AF ≥ threshold (0.6) | **ALT 염기로 확정** | 지배적인 변이 allele |
| AF ≤ 1 − threshold (0.4) | **REF 염기로 확정** | 지배적인 reference allele |
| 0.4 < AF < 0.6 | **IUPAC 유지** | 혼합감염 의심 |

> threshold는 `ivar_threshold` 파라미터(기본 0.6)를 공유합니다.

#### AF 정보 파일 (`*.af.tsv`)

`results/translation/per_sample/{barcode}.{GENE}.af.tsv` 에 저장됩니다.

| 컬럼 | 설명 |
|---|---|
| `sample`, `gene` | 샘플 / 유전자 |
| `chrom`, `gen_pos` | 게놈 위치 (1-based) |
| `cds_pos`, `codon`, `frame` | CDS 내 위치 / 코돈 번호 / 코돈 내 위치(0~2) |
| `iupac` | 원래 IUPAC 코드 (Y, R, K 등) |
| `ref`, `alt` | VCF 상의 REF / ALT 염기 |
| `ref_depth`, `alt_depth`, `total_depth` | 각 allele 읽기 수 |
| `af` | ALT allele frequency (alt / total) |
| `resolved` | 최종 결정된 염기 |
| `note` | `alt_call` / `ref_call` / `ambiguous_mixed` |

**`ambiguous_mixed` 행이 많은 샘플은 혼합감염(mixed infection)을 의심**할 수 있으며,
`alt_call` / `ref_call`로 해소된 경우 번역 결과에서 X 대신 실제 아미노산이 출력됩니다.

---

## ⚙️ 설정 파일 우선순위

파이프라인은 두 개의 설정 파일을 사용합니다:

| 파일 | 자동 적용 | 용도 |
|---|---|---|
| `nextflow.config` | ✅ **자동** (항상 적용) | 모든 파라미터 기본값, 실행 프로파일 |
| `configs/params.yaml` | ❌ **수동** (`-params-file` 옵션 필요) | 파라미터 참고·설명 문서 |

`configs/params.yaml`을 실제 실행에 반영하려면:

```bash
nextflow run main.nf -profile standard -params-file configs/params.yaml
```

> 일반적으로는 `nextflow.config`에서 직접 수정하는 것이 가장 간단합니다.

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
